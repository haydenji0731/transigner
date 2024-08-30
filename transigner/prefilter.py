#!/usr/bin/env python

from transigner.commons import *
from transigner.model import *

def load_txome_info(fh):
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading transcriptome information")
    tnames = fh.header.references
    tlens = fh.header.lengths
    assert len(tnames) == len(tlens)
    ti_map = dict()
    for ti in range(len(tnames)):
        ti_map[tnames[ti]] = ti # tname - ti mapping
    return ti_map, tlens

def load_alignments(fn):
    acore_mat = []
    ctr = -1
    qi_map = dict() # massive performance improvement
    pri_tis = []
    unmapped = set()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        ti_map, tlens = load_txome_info(fh)
        ti_size = len(ti_map)
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading alignment records")
        for brec in tqdm(fh):
            if brec.is_unmapped:
                unmapped.add(brec.query_name) # excluded from qi_map
                continue
            elif brec.is_supplementary:
                continue
            else:
                qname = brec.query_name
                if qname not in qi_map:
                    ctr += 1
                    qi_map[qname] = ctr
                    acore_mat.append(dict())
                    # acore_mat.append([None] * ti_size)
                    pri_tis.append(-1)
                    qi = ctr
                else:
                    qi = qi_map[qname]
                ti = ti_map[brec.reference_name]
                curr_acore = acore(rst=brec.reference_start, \
                                   ren=brec.reference_end, \
                                   ascore=int(brec.get_tag("ms")), \
                                   pscore=-1)
                if ti in acore_mat[qi]:
                    if curr_acore.ascore > acore_mat[qi][ti].ascore:
                        acore_mat[qi][ti] = curr_acore
                else:
                    acore_mat[qi][ti] = curr_acore
                
                if not brec.is_secondary:
                    assert pri_tis[qi] == -1
                    pri_tis[qi] = ti
    
    return acore_mat, pri_tis, qi_map, ti_map, tlens, unmapped

def calc_pscores(qi_size, tlens, acore_mat, resume=False):
    if resume:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} resuming a previous run")
        with open('tsegs_mat.pkl', 'rb') as f:
            tsegs_mat = pickle.load(f)
        with open('n_lst.pkl', 'rb') as f:
            n_lst = pickle.load(f)
        with open('pspb_mat.pkl', 'rb') as f:
            log_pspb_mat = pickle.load(f)
    else:
        tsegs_mat, n_lst = get_scovs(qi_size, tlens, acore_mat)
        log_pspb_mat = calc_log_pspbs(tsegs_mat, tlens, n_lst)

    for qi in tqdm(range(qi_size)):
        for ti in acore_mat[qi]:
            assert ti in log_pspb_mat
            curr_acore = acore_mat[qi][ti]
            pscore = get_pscore(curr_acore, ti, log_pspb_mat)
            new_acore = curr_acore._replace(pscore=pscore)
            acore_mat[qi][ti] = new_acore

def build_cmpt_mat(acore_mat, qi_size, pri_tis):
    cmpt_mat = [dict() for _ in range(qi_size)]
    for qi in tqdm(range(qi_size)):
        pri_ti = pri_tis[qi]
        if pri_ti == -1:
            print(datetime.now(), f"{RED}ERROR{RESET} no primary alignment for a read")
            sys.exit(-1)
        pri_acore = acore_mat[qi][pri_ti]
        cmpt_mat[qi][pri_ti] = 1.0 * pri_acore.pscore
        for ti in acore_mat[qi]:
            curr_acore = acore_mat[qi][ti]
            sigma_a = math.exp((curr_acore.ascore - pri_acore.ascore) / 10)
            sigma_p = curr_acore.pscore
            sigma_comb = sigma_a * sigma_p
            if cmpt_mat[qi][ti] is not None:
                cmpt_mat[qi][ti] = max(cmpt_mat[qi][ti], sigma_comb)
            else:
                cmpt_mat[qi][ti] = sigma_comb
    return cmpt_mat

def write_cmpt_mat(cmpt_mat, qi_map, ti_map, qi_size, out_dir):
    cmpt_lns = []
    score_lns = []
    qnames = sorted(qi_map, key=lambda k: qi_map[k])
    tnames = sorted(ti_map, key=lambda k: ti_map[k])
    for qi in tqdm(range(qi_size)):
        qname = qnames[qi]
        tnames_qi = []
        scores_qi = []
        for ti in cmpt_mat[qi]:
            tname = tnames[ti]
            tnames_qi.append(tname)
            scores_qi.append(f'({tname},{cmpt_mat[qi][ti]})')
        cmpt_lns.append(f'{qname}' + '\t'.join(tnames_qi))
        score_lns.append(f'{qname}' + '\t'.join(scores_qi))
    
    cmpt_fn = os.path.join(out_dir, 'cmpt_mat.tsv')
    with open(cmpt_fn, 'w') as f:
        f.write('\n'.join(cmpt_lns) + '\n')
    
    score_fn = os.path.join(out_dir, 'scores.tsv')
    with open(score_fn, 'w') as f:
        f.write('\n'.join(score_lns) + '\n')

def main(args):
    cmd_fn = os.path.join(args.out_dir, "pref_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)
    
    acore_mat, pri_tis, qi_map, ti_map, tlens, unmapped = load_alignments(args.aln)
    qi_size = len(qi_map)
    ti_size = len(ti_map)

    unmapped_lns = []
    for qname in unmapped:
        unmapped_lns.append(f'{qname}')
    unm_fn = os.path.join(args.out_dir, "unmapped.txt")
    with open(unm_fn, 'w') as f:
        f.write('\n'.join(unmapped_lns) + "\n")

    # saves ti dict
    ti_ifn = os.path.join(args.out_dir, "ti.pkl")
    with open(ti_ifn, 'wb') as f:
        pickle.dump(ti_map, f)
    
    if args.tsv:
        qi_tfn = os.path.join(args.out_dir, 'qi_map.tsv')
        to_tsv(qi_map, qi_tfn)
        ti_tfn = os.path.join(args.out_dir, 'ti_map.tsv')
        to_tsv(ti_map, ti_tfn)
    
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} computing position-specific scores")
    calc_pscores(qi_size, tlens, acore_mat, args.resume)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} building compatibility matrix")
    cmpt_mat = build_cmpt_mat(acore_mat, qi_size, pri_tis)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out")
    write_cmpt_mat(cmpt_mat, qi_map, ti_map, qi_size, args.out_dir)

if __name__ == "__main__":
    main()