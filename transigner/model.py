from transigner.commons import *

def generate_tsegs(tsegs_mat, ti, tlen, n_seg):
    seg_len = math.floor(tlen / n_seg)
    adj_n_seg = math.ceil(tlen / seg_len)
    tsegs = [None] * adj_n_seg
    for si in range(adj_n_seg):
        sst = si * seg_len
        sen = min((si + 1) * seg_len, tlen)
        slen = sen - sst
        tseg_si = tseg(parent=ti, sst=sst, sen=sen, slen=slen, scov=0, si=si)
        tsegs[si] = tseg_si
    tsegs_mat[ti] = tsegs
    return adj_n_seg

# acore_tbl; <k,v> = <ti, list(acore)>
def get_acores_by_ti(acore_mat, qi_size):
    acores_by_ti = dict()
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} aggregating alignment records by transcript")
    for qi in tqdm(range(qi_size)):
        for ti in acore_mat[qi]:
            if ti not in acores_by_ti:
                acores_by_ti[ti] = [acore_mat[qi][ti]]
            else:
                acores_by_ti[ti].append(acore_mat[qi][ti])
    return acores_by_ti


def get_pbase_cov(ti, tlen, acores_by_ti):
    pbase_cov = [0] * tlen
    acores_ti = acores_by_ti[ti]
    for acore in acores_ti:
        for i in range(acore.rst, acore.ren, 1):
            pbase_cov[i] += 1
    return pbase_cov

# acore_mat = list of dict()
def get_scovs(qi_size, tlens, acore_mat, n_seg=10):
    tsegs_mat = dict()
    n_lst = dict()
    acores_by_ti = get_acores_by_ti(acore_mat, qi_size)
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} computing per-segment coverages")
    for ti in tqdm(acores_by_ti):
        tlen = tlens[ti]
        adj_n_seg = generate_tsegs(tsegs_mat, ti, tlen, n_seg) # either n_seg or n_seg + 1
        # sanity check
        assert tsegs_mat[ti] is not None
        assert len(tsegs_mat[ti]) == adj_n_seg
        pbase_cov = get_pbase_cov(ti, tlen, acores_by_ti)
        n_lst[ti] = 0
        scov = 0
        for si in range(adj_n_seg):
            tseg_si = tsegs_mat[ti][si]
            try:
                scov = sum(pbase_cov[tseg_si.sst:tseg_si.sen]) / tseg_si.slen
            except:
                print(tseg_si)
                sys.exit(-1)
            new_tseg_si = tseg_si._replace(scov=scov)
            tsegs_mat[ti][si] = new_tseg_si
            n_lst[ti] += scov
    with open("tsegs_mat.pkl", 'wb') as f:
        pickle.dump(tsegs_mat, f)
    
    with open("n_lst.pkl", 'wb') as f:
        pickle.dump(n_lst, f)
    return tsegs_mat, n_lst

def calc_log_sspb(n, k_i):
    p_i = k_i / n
    gam_n = gammaln(n + 1)
    gam_k_i = gammaln(k_i + 1)
    gam_n_k_i = gammaln(n - k_i + 1)
    log_sspb = (gam_n - gam_k_i - gam_n_k_i)
    if n == k_i:
        log_sspb += k_i * np.log(p_i)
    else:
        log_sspb += (k_i * np.log(p_i)) + ((n - k_i) * np.log(1 - p_i))
    return log_sspb

# computes position-specific probabilities; p of covering a base
def calc_log_pspbs(tsegs_mat, tlens, n_lst):
    log_pspb_mat = dict()
    for ti in tqdm(tsegs_mat):
        tlen = tlens[ti]
        tsegs = tsegs_mat[ti]
        n_seg = len(tsegs) # transcript-specific number of segments
        n_ti = n_lst[ti]
        log_pspb_ti = [None] * tlen # position-specific probability
        log_sspbs = [None] * n_seg
        for si in range(n_seg):
            tseg_si = tsegs[si]
            k_si = tseg_si.scov
            if k_si == 0:
                log_sspbs[si] = float('nan')
                continue
            log_sspbs[si] = calc_log_sspb(n_ti, k_si) #segment-specific probability
        # sanity check
        assert None not in log_sspbs

        min_log_sspb = abs(min([x for x in log_sspbs if not math.isnan(x)])) + 1e-3 # does the padding matter?
        log_sspbs_shift = [min_log_sspb + x if not math.isnan(x) else x for x in log_sspbs]
        log_sspbs_shift = [0 if math.isnan(x) else x for x in log_sspbs_shift]

        # normalize sspbs
        log_sspbs_sum = sum(log_sspbs_shift)
        norm_log_sspbs = [p / log_sspbs_sum for p in log_sspbs_shift]
        for si in range(n_seg):
            tseg_si = tsegs[si]
            log_pspb = norm_log_sspbs[si] / tseg_si.slen
            for pi in range(tseg_si.sst, tseg_si.sen, 1):
                log_pspb_ti[pi] = log_pspb
        # sanity check
        assert None not in log_pspb_ti
        log_pspb_mat[ti] = log_pspb_ti

    # for debugging
    with open("pspb_mat.pkl", 'wb') as f:
        pickle.dump(log_pspb_mat, f)
    return log_pspb_mat

def get_pscore(acore, ti, pspb_mat):
    pscore = 0
    for pi in range(acore.rst, acore.ren, 1):
        pscore += pspb_mat[ti][pi]
    return pscore