from transigner.utils import *
import pysam
import numpy as np
import pandas as pd
import math
from collections import namedtuple

aObj = namedtuple('aObj', ['qname', 'tname', 'start', 'end', 'score'])

def load_bam(fn):
    amat = dict()
    unmapped = set()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        print(tmessage("retrieving reference info", Mtype.PROG))
        tlens = fh.header.lengths
        tnames = fh.header.references
        tmap = dict()
        for i, tname in enumerate(tnames):
            tmap[tname] = i
        print(tmessage("loading bam records", Mtype.PROG))
        for brec in fh:
            if brec.is_unmapped: unmapped.add(brec.query_name); continue
            elif brec.is_supplementary: continue
            else:
                qname = brec.query_name
                tname = brec.reference_name
                if qname not in amat:
                    amat[qname] = dict()
                ti = tmap[tname]
                aobj = aObj(qname=qname, tname=brec.reference_name, \
                            start=brec.reference_start, end=brec.reference_end, \
                            score=int(brec.get_tag("AS")))
                if ti in amat[qname]:
                    amat[qname][ti] = aobj if aobj.score > amat[qname][ti].score else amat[qname][ti]
                else:
                    amat[qname][ti] = aobj
    return np.array(tlens), tmap, amat, unmapped

def save_tmap(tmap, fn):
    out_lns = []
    for tname in tmap:
        out_lns.append(f'{tname},{tmap[tname]}')
    with open(fn, 'w') as fh:
        fh.write('\n'.join(out_lns) + '\n')

def load_estimates(fn) -> dict:
    df = pd.read_csv(fn)
    est_tbl = dict()
    for i, row in df.iterrows():
        tid = row['transcript_id']
        cov = float(row['cov'])
        fpkm = float(row['fpkm'])
        tpm = float(row['tpm'])
        est_tbl[tid] = (cov, fpkm, tpm)
    return est_tbl
        
def calc_psw(amat, tlens):
    tsize = len(tlens)
    per_tx_cov = np.zeros(tsize)
    per_base_cov = [np.zeros(l) for l in tlens]
    for qname in amat:
        for ti in amat[qname]:
            aobj = amat[qname][ti]
            per_tx_cov[ti] += aobj.end - aobj.start
            per_base_cov[ti][aobj.start:aobj.end] += 1
    per_tx_cov /= tlens
    psw = dict()
    for ti in range(tsize):
        if per_tx_cov[ti] == 0: continue
        tlen = tlens[ti]
        delt = per_tx_cov[ti] - per_base_cov[ti]
        delt += np.abs(np.min(delt))
        z = np.sum(delt)
        if z == 0: psw[ti] = np.ones(tlen) / tlen; continue
        psw[ti] = delt / z
    cum_psw = dict()
    for ti in psw:
        cum_psw[ti] = np.cumsum(psw[ti])
    return per_tx_cov, per_base_cov, psw, cum_psw

def calc_psw_pre(amat, est, tlens, tnames):
    tsize = len(tlens)
    per_base_cov = [np.zeros(l) for l in tlens]
    for qname in amat:
        for ti in amat[qname]:
            aobj = amat[qname][ti]
            per_base_cov[ti][aobj.start:aobj.end] += 1
    psw = dict()
    for ti in range(tsize):
        tname = tnames[ti]
        assert tname in est # sanity check
        tlen = tlens[ti]
        delt = est[tname][0] - per_base_cov[ti]
        delt += np.abs(np.min(delt))
        z = np.sum(delt)
        if z == 0: psw[ti] = np.ones(tlen) / tlen; continue
        psw[ti] = delt / z
    cum_psw = dict()
    for ti in psw:
        cum_psw[ti] = np.cumsum(psw[ti])
    return per_base_cov, psw, cum_psw

def build_cmat(amat, cum_psw, opts_d):
    cmat = {x:dict() for x in amat}
    for qname in amat:
        assert qname in cmat # sanity check
        max_score = max([x.score for x in list(amat[qname].values())])
        z = 0
        temp_1 = []
        temp_2 = []
        for ti in amat[qname]:
            aobj = amat[qname][ti]
            score_delt = max_score - aobj.score
            if opts_d['score_mdl'] == 'e':
                sigma_a = math.exp((-1 * score_delt) / opts_d['score_K']) # TODO: check if this has an impact
            else:
                sigma_a = opts_d['score_A'] / score_delt ** opts_d['score_K']
            if opts_d['use_psw']:
                sigma_p = cum_psw[ti][aobj.end - 1] - cum_psw[ti][aobj.start]
            else:
                sigma_p = 1.0
            sigma = sigma_a * sigma_p
            z += sigma
            temp_1.append((ti, sigma_a))
            temp_2.append((ti, sigma))
        if z > 0:
            cmat[qname] = dict(temp_2)
        else:
            cmat[qname] = dict(temp_1)
        # cmat[qname][ti] = sigma_a * sigma_p
    return cmat

def load_opts(args):
    opts_d = dict()
    opts_d['score_mdl'] = args.dp_score_model
    opts_d['score_A'] = args.dsm_opts[0]
    opts_d['score_K'] = args.dsm_opts[1]
    opts_d['use_filt'] = args.use_filter
    opts_d['filt_fp'] = args.filt_opts[0]
    opts_d['filt_tp'] = args.filt_opts[1]
    opts_d['filt_tcov'] = args.filt_opts[2]
    opts_d['use_psw'] = args.use_psw
    return opts_d

def join_and_write(lst, fn):
    with open(fn, 'w') as fh: fh.write('\n'.join(lst) + "\n")

# TODO: debug
def write_psw_data(tsize, per_tx_cov, per_base_cov, psw, cum_psw, tnames, out_dir):
    fn = os.path.join(out_dir, 'per_transcript_cov.csv')
    out_lns = [f'{tnames[i],{per_tx_cov[i]}}' for i in range(tsize)]
    join_and_write(out_lns, fn)
    out_lns_pbase = []
    out_lns_psw = []
    out_lns_cpsw = []
    for i in range(tsize):
        tname = tnames[i]
        out_lns_pbase.append(f'{tname},{np.array2string(per_base_cov[i], separator=' ')}')
        out_lns_psw.append(f'{tname},{np.array2string(psw[i], separator=' ')}')
        out_lns_cpsw.append(f'{tname},{np.array2string(cum_psw[i], separator=' ')}')
    fn = os.path.join(out_dir, 'per_base_cov.csv')
    join_and_write(out_lns_pbase, fn)
    fn = os.path.join(out_dir, 'psw.csv')
    join_and_write(out_lns_psw, fn)
    fn = os.path.join(out_dir, 'cum_psw.csv')
    join_and_write(out_lns_cpsw, fn)

def write_cmat(cmat, tnames, fn):
    out_lns = []
    for qname in cmat:
        for ti in cmat[qname]:
            out_lns.append(f'{qname},{tnames[ti]},{cmat[qname][ti]}')
    with open(fn, 'w') as fh:
        fh.write(f'#{len(cmat)};{len(tnames)}\n')
        fh.write('\n'.join(out_lns) + "\n")
    
def write_tp_scores(cmat, tnames, fn):
    out_lns = []
    for qname in cmat:
        z = sum(cmat[qname].values())
        if z == 0: continue
        for ti in cmat[qname]:
            tname = tnames[ti]
            if '_'.join(qname.split('_')[0:2]) == tname.split('.')[0]:
                out_lns.append(f'{qname},{tname},{cmat[qname][ti] / z}')
    with open(fn, 'w') as fh:
        fh.write(f'read_id,transcript_id,score\n')
        fh.write('\n'.join(out_lns) + "\n")

def main(args) -> None:
    check_dir(args.out_dir)
    param_fn = os.path.join(args.out_dir, "pre_params.json")
    store_params(args, param_fn)
    opts_d = load_opts(args)

    tlens, tmap, amat, unmapped = load_bam(args.in_file)
    unmapped_fn = os.path.join(args.out_dir, "unmapped.txt")
    with open(unmapped_fn, 'w') as f:
        for x in unmapped: f.write(f'{x}\n')

    tnames = sorted(tmap, key=lambda k: tmap[k])

    tmap_fn = os.path.join(args.out_dir, "tmap.csv")
    save_tmap(tmap, tmap_fn)
    
    use_psw = args.use_psw if not args.spiked else False
    opts_d['use_psw'] = use_psw
    
    if use_psw:
        if not args.estimates:
            print(tmessage("calculating position specific weights", Mtype.PROG))
            per_tx_cov, per_base_cov, psw, cum_psw = calc_psw(amat, tlens)
        else:
            print(tmessage("loading coverage estimates", Mtype.PROG))
            est_tbl = load_estimates(args.estimates)
            per_base_cov, psw, cum_psw = calc_psw_pre(amat, est_tbl, tlens, tnames)
    else:
        cum_psw = None

    print(tmessage("building compatibility matrix", Mtype.PROG))
    cmat = build_cmat(amat, cum_psw, opts_d)
    
    score_fn = os.path.join(args.out_dir, "scores.csv")
    write_cmat(cmat, tnames, score_fn)

    if args.dev:
        score_fn = os.path.join(args.out_dir, "tp_scores.csv")
        write_tp_scores(cmat, tnames, score_fn)
