#!/usr/bin/env python

import pysam
from tqdm import tqdm
import os
import pickle
from collections import namedtuple
import json
from datetime import datetime
from transigner.commons import RED, GREEN, RESET
import numpy as np
import math
import sys

acore = namedtuple('acore', ['rst', 'ren', 'score'])

def load_aln(fn):
    aln_mat = []
    ctr = -1
    qi_map = dict()
    pri_lst = []
    unmapped = set()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} building target transcriptome index")
        tlens = fh.header.lengths
        tnames = fh.header.references
        ti_map = dict()
        for i, tname in enumerate(tnames):
            ti_map[tname] = i
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading alignment records")
        for brec in fh:
            if brec.is_unmapped:
                unmapped.add(brec.query_name)
                continue
            elif brec.is_supplementary:
                continue
            else:
                qname = brec.query_name
                if qname not in qi_map:
                    ctr += 1
                    qi_map[qname] = ctr
                    aln_mat.append(dict())
                    pri_lst.append(-1)
                    qi = ctr
                else:
                    qi = qi_map[qname]
                ti = ti_map[brec.reference_name]
                # NOTE: reference_end is 1 past the last aligned residue
                curr_acore = acore(rst=brec.reference_start, \
                                    ren=brec.reference_end, \
                                    score=int(brec.get_tag("AS"))) 
                if ti in aln_mat[qi]:
                    if curr_acore.score > aln_mat[qi][ti].score:
                        aln_mat[qi][ti] = curr_acore
                else:
                    aln_mat[qi][ti] = curr_acore
                # NOTE: primary alignment is assigned based on ms
                # src: https://github.com/lh3/minimap2/issues/412
                if not brec.is_secondary:
                    assert pri_lst[qi] == -1
                    pri_lst[qi] = ti
    return tlens, tnames, ti_map, aln_mat, pri_lst, qi_map, unmapped

def calc_psw(aln_mat, qi_size, tlens, tnames):
    nb_tbl = dict()
    pb_cov_tbl = dict()
    for qi in range(qi_size):
        for ti in aln_mat[qi]:
            acore_curr = aln_mat[qi][ti]
            if ti not in nb_tbl:
                nb_tbl[ti] = acore_curr.ren - acore_curr.rst
            else:
                nb_tbl[ti] += acore_curr.ren - acore_curr.rst
            if ti not in pb_cov_tbl:
                pb_cov_tbl[ti] = np.zeros(tlens[ti])
            pb_cov_tbl[ti][acore_curr.rst:acore_curr.ren] += 1

    pt_cov_tbl = dict()
    psw_tbl = dict()
    for ti in tqdm(nb_tbl):
        pt_cov_tbl[ti] = nb_tbl[ti] / tlens[ti]
        delt = pt_cov_tbl[ti] - pb_cov_tbl[ti]
        delt += np.abs(np.min(delt))
        z = np.sum(delt) # TODO: think about this part a little harder
        if z == 0:
            psw_tbl[ti] = np.ones(tlens[ti])
            psw_tbl[ti] /= tlens[ti]
            continue
        psw_tbl[ti] = delt / z
    return nb_tbl, pt_cov_tbl, pb_cov_tbl, psw_tbl
    
def build_cmpt_mat(aln_mat, qi_size, pri_lst, tlens, cum_psw_tbl, \
                   filter_args, use_psw=False, use_filter=False):
    if use_filter:
        if not filter_args:
            raise Exception(f"Provide filter args")
        fp_thres, tp_thres = filter_args

    cmpt_mat = [dict() for _ in range(qi_size)]
    for qi in tqdm(range(qi_size)):
        max_score = max([x.score for x in list(aln_mat[qi].values())])
        for ti in aln_mat[qi]:
            curr_acore = aln_mat[qi][ti]
            if curr_acore.score == max_score:
                sigma_a = 1.0
            else:
                x = max_score - curr_acore.score
                sigma_a = 1 / x
            # sigma_a = math.exp((curr_acore.score - max_score))
            if use_psw:
                sigma_p = cum_psw_tbl[ti][curr_acore.ren - 1] - cum_psw_tbl[ti][curr_acore.rst]
            else:
                sigma_p = 1.0
            cmpt_mat[qi][ti] = sigma_a * sigma_p

    # for qi in tqdm(range(qi_size)):
    #     pri_ti = pri_lst[qi]
    #     if pri_ti == -1:
    #         raise Exception(f"No primary alignment for read index {qi}")
        
    #     pri_rst, pri_ren, pri_score = aln_mat[qi][pri_ti]
    #     if use_psw:
    #         # if pri_rst == 0:
    #         #     sigma_p = cum_psw_tbl[pri_ti][pri_ren - 1]
    #         # else:
    #         sigma_p = cum_psw_tbl[pri_ti][pri_ren - 1] - cum_psw_tbl[pri_ti][pri_rst]
    #     else:
    #         sigma_p = 1.0
    #     cmpt_mat[qi][pri_ti] = 1.0 * sigma_p
        
    #     for ti in aln_mat[qi]:
    #         if ti == pri_ti:
    #             continue
    #         acore = aln_mat[qi][ti]

    #         if use_filter:
    #             fp_dist = pri_rst - acore.rst
    #             tlen = tlens[ti]
    #             tp_dist = (tlens[pri_ti] - pri_ren) - (tlen - acore.ren)
    #             if fp_dist < fp_thres or tp_dist < tp_thres:
    #                 continue # filter
    
    #         # x = (acore.score - pri_score) / 5
    #         # for safety
    #         # if x > 709:
    #         #     x = 709
    #         # elif x < -709:
    #         #     x = -709
    #         # sigma_a = math.exp(x)
    #         sigma_a = acore.score / pri_score

    #         if use_psw:
    #             # if acore.rst == 0:
    #             #     sigma_p = cum_psw_tbl[ti][acore.ren - 1]
    #             # else:
    #             sigma_p = cum_psw_tbl[ti][acore.ren - 1] - cum_psw_tbl[ti][acore.rst]
    #         else:
    #             sigma_p = 1.0
    #         # cmpt_mat[qi][ti] = sigma_a * sigma_p

    #         cmpt_score = sigma_a * sigma_p
    #         if ti in cmpt_mat[qi]:
    #             cmpt_mat[qi][ti] = max(cmpt_mat[qi][ti], cmpt_score)
    #         else:
    #             cmpt_mat[qi][ti] = cmpt_score
    return cmpt_mat

def main(args) -> None:
    cmd_fn = os.path.join(args.out_dir, "pref_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    tlens, tnames, ti_map, aln_mat, pri_lst, qi_map, unmapped = load_aln(args.aln)

    unm_fn = os.path.join(args.out_dir, "unmapped.txt")
    unm_fh = open(unm_fn, 'w')
    for qname in unmapped:
        unm_fh.write(qname + "\n")
    unm_fh.close()

    ti_map_fn = os.path.join(args.out_dir, "ti.pkl")
    with open(ti_map_fn, "wb") as f:
        pickle.dump(ti_map, f)

    if args.psw:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} calculating position specific weights")
        nb_tbl, pt_cov_tbl, pb_cov_tbl, psw_tbl = calc_psw(aln_mat, len(qi_map), tlens, tnames)
        cum_psw_tbl = dict()
        for ti in tqdm(psw_tbl):
            cum_psw_tbl[ti] = np.cumsum(psw_tbl[ti])
        # print(cum_psw_tbl[5][2239] - cum_psw_tbl[5][36])
        # print(cum_psw_tbl[5][2240] - cum_psw_tbl[5][36])
        # print(cum_psw_tbl[5][2240])
        # print(cum_psw_tbl[5][2239])
        # print(cum_psw_tbl[5][36])
        # print(len(cum_psw_tbl[5]))
        # sys.exit(-1)
    else:
        cum_psw_tbl = None

        
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} building a compatibility matrix")
    cmpt_mat = build_cmpt_mat(aln_mat, len(qi_map), pri_lst, tlens, \
                            cum_psw_tbl, args.filt_opts, args.psw, \
                            args.filter)
    
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out")
    cmpt_fn = os.path.join(args.out_dir, "cmpt_mat.csv")
    cmpt_fh = open(cmpt_fn, 'w')
    cmpt_fh.write("#" + str(len(cmpt_mat)) + "\n")
    score_fn = os.path.join(args.out_dir, "scores.csv")
    score_fh = open(score_fn, 'w')
    score_fh.write("#" + str(len(cmpt_mat)) + "\n")
    sur_fn = os.path.join(args.out_dir, "unassigned.txt")
    sur_fh = open(sur_fn, 'w')
    qnames = sorted(qi_map, key=lambda k: qi_map[k])
    tnames = sorted(ti_map, key=lambda k: ti_map[k])
    # TODO: make this write faster
    for qi in range(len(qi_map)):
        qname = qnames[qi]
        if len(cmpt_mat[qi]) == 0:
            sur_fh.write(qname + "\n")
            continue
        for ti in cmpt_mat[qi]:
            tname = tnames[ti]
            cmpt_fh.write(f'{qname},{tname}\n')
            score_fh.write(f'{qname},{tname},{cmpt_mat[qi][ti]}\n')
    sur_fh.close()
    cmpt_fh.close()
    score_fh.close()

if __name__ == "__main__":
    main()