#!/usr/bin/env python

import pysam
from tqdm import tqdm
import pyfastx
import os
import pickle
from collections import namedtuple
import sys
import json
from datetime import datetime
from transigner.commons import RED, GREEN, RESET
from joblib import Parallel, delayed

# model fitting imports
import numpy as np
from scipy.optimize import minimize
from scipy.stats import gaussian_kde, gamma
import scipy.stats as stats
from scipy.signal import find_peaks
from diptest import diptest

acore = namedtuple('acore', ['rst', 'ren', 'score'])
gpdf = namedtuple('gpdf', ['alpha', 'beta'])

def load_aln(fn, ti_map, model_pos=False):
    aln_mat = []
    ctr = -1
    qi_map = dict()
    pri_lst = []
    unmapped = set()
    rpos_tbl = None
    if model_pos:
        rpos_tbl = dict()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in tqdm(fh):
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
                curr_acore = acore(rst=brec.reference_start, \
                                    ren=brec.reference_end, \
                                    score=int(brec.get_tag("ms")))
                if ti in aln_mat[qi]:
                    if curr_acore.score > aln_mat[qi][ti].score:
                        aln_mat[qi][ti] = curr_acore
                else:
                    aln_mat[qi][ti] = curr_acore

                if not brec.is_secondary:
                    assert pri_lst[qi] == -1
                    pri_lst[qi] = ti
                # for now, only modeling start pos in dRNA reads
                if model_pos:
                    if ti not in rpos_tbl:
                        rpos_tbl[ti] = [brec.reference_start]
                    else:
                        rpos_tbl[ti].append(brec.reference_start)
    fh = open("qi_map.tsv", 'w')
    for qname in qi_map:
        fh.write(qname + "\t" + str(qi_map[qname]) + "\n")
    fh.close()
    return aln_mat, pri_lst, qi_map, unmapped, rpos_tbl

def build_ti_map(fn):
    fa = pyfastx.Fasta(fn)
    ti_map = dict()
    tlen_lst = []
    ctr = 0
    for obj in fa:
        ti_map[obj.name] = ctr
        tlen_lst.append(len(obj.seq))
        ctr += 1
    fh = open("ti_map.tsv", 'w')
    for tname in ti_map:
        fh.write(tname + "\t" + str(ti_map[tname]) + "\n")
    fh.close()
    return ti_map, tlen_lst

def build_cmpt_mat(aln_mat, qi_map, pri_lst, tlen_lst, \
                filter, fp_thres, tp_thres, surrender, \
                tcov_thres, pos_dist, model=False):
    qi_size = len(qi_map)
    cmpt_mat = [dict() for _ in range(qi_size)]
    sur_ctr = 0
    out_fh = open("./pos_scores.tsv", 'w')
    print("wa")
    for qi in tqdm(range(qi_size)):
        pri_ti = pri_lst[qi]
        if pri_ti == -1:
            print("FATAL ERROR: no primary aln for a read")
            sys.exit(-1)
        pri_rst, pri_ren, pri_score = aln_mat[qi][pri_ti]
        pri_tlen = tlen_lst[pri_ti]
        tcov = (pri_ren - pri_rst) / tlen_lst[pri_ti]
        if surrender:
            tcov = (pri_ren - pri_rst) / tlen_lst[pri_ti]
            # surrender reads with poor tcov
            if tcov < tcov_thres:
                continue
        if model:
            if pri_ti in pos_dist:
                pos_score = gamma.pdf(pri_rst, a=pos_dist[pri_ti].alpha, scale=pos_dist[pri_ti].beta)
                if np.isinf(pos_score):
                    pos_score = 1.0
            else:
                pos_score = None
            cmpt_mat[qi][pri_ti] = 1.0
        else:
            pos_score = Nones
            cmpt_mat[qi][pri_ti] = 1.0
        out_fh.write(f'{qi}\t{pri_ti}\t1.0\t{pos_score}\t{tcov}\n')

        for ti in aln_mat[qi]:
            if ti == pri_ti:
                continue
            acore = aln_mat[qi][ti]
            tlen = tlen_lst[ti]
            tcov = (acore.ren - acore.rst) / tlen
            if filter: # no filter performs better?
                fp_dist = pri_rst - acore.rst
                if tp_thres == -1:
                    if fp_dist < fp_thres:
                        continue
                tlen = tlen_lst[ti]
                tp_dist = (pri_tlen - pri_ren) - (tlen - acore.ren)
                if fp_dist < fp_thres or tp_dist < tp_thres:
                    continue
            
            if model:
                if ti in pos_dist:
                    pos_score = gamma.pdf(acore.rst, a=pos_dist[ti].alpha, scale=pos_dist[ti].beta)
                    if np.isinf(pos_score):
                        pos_score = 1.0
                else:
                    pos_score = None
            else:
                pos_score = None
            cmpt_score = acore.score / pri_score
            # cmpt_score = acore.score / pri_score + 0.05 * pos_score + 0.1 * tcov
            out_fh.write(f'{qi}\t{ti}\t{cmpt_score}\t{pos_score}\t{tcov}\n')
            if ti in cmpt_mat[qi]:
                cmpt_mat[qi][ti] = max(cmpt_mat[qi][ti], cmpt_score)
            else:
                cmpt_mat[qi][ti] = cmpt_score
        if len(cmpt_mat[qi]) == 0:
            sur_ctr += 1
    out_fh.close()
    return cmpt_mat, sur_ctr

def count_peaks(data, bw=0.5, h=0.01):
    kde = gaussian_kde(data, bw_method=bw)
    x = np.linspace(min(data), max(data), 1000)
    y = kde.evaluate(x)
    peaks, _ = find_peaks(y, height=h)
    return len(peaks)

def gamma_mixture_ll(params, x):
    a1, b1, a2, b2, w1, w2 = params
    if not (0 < w1 < 1 and 0 < w2 < 1 and np.isclose(w1 + w2, 1)):
        return np.inf
    if a1 <= 0 or b1 <= 0 or a2 <= 0 or b2 <= 0:
        return np.inf
    pdf1 = gamma.pdf(x, a1, scale=b1)
    pdf2 = gamma.pdf(x, a2, scale=b2)
    mixture_pdf = w1 * pdf1 + w2 * pdf2
    return -np.sum(np.log(mixture_pdf + 1e-10))

def fit_gamma(data):
    bin_w, n_bins = freedman_diaconis(data)
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2
    initial_guess = [2.0, 2.0, 0.5, 3.0, 0.5, 0.5]
    bounds = [
    (0.1, None),
    (0.1, None),
    (0, 1),
    (0.1, None),
    (0.1, None),
    (0, 1)
    ]
    constraints = [{'type': 'eq', 'fun': lambda params: np.sum(params[4:]) - 1}]
    result = minimize(gamma_mixture_ll, initial_guess, args=(x,), bounds=bounds, \
    constraints=constraints, method='SLSQP')
    return result

def get_tname(ti_map, ti):
    for tname in ti_map:
        if ti_map[tname] == ti:
            return tname
    return None

def freedman_diaconis(data):
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_w = 2 * iqr * len(data) ** (-1/3)
    n_bins = int(np.ceil(data.max() - data.min()) / bin_w)
    return bin_w, n_bins

def main(args):
    cmd_fn = os.path.join(args.out_dir, "pref_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    # TODO: speed this up by just looking at the bam hdr
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading target transcriptome")
    ti_map, tlen_lst = build_ti_map(args.target)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading query-to-target alignments")
    aln_mat, pri_lst, qi_map, unmapped, rpos_tbl = load_aln(args.aln, ti_map, args.model) # rpos_tbl is None if not args.model
    # rpos_tbl stores reference start positions; currently only modeling aln start pos

    if args.model:
        dpeak_fn = os.path.join(args.out_dir, "double_peaks.txt")
        dpeak_fh = open(dpeak_fn, 'w')
        pos_dist = dict()
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} checking multimodality in alignment positions")
        for ti in tqdm(rpos_tbl):
            rpos_arr = np.array(rpos_tbl[ti])
            tname = get_tname(ti_map, ti)
            if len(rpos_tbl[ti]) < 1000:
                continue
            elif len(rpos_tbl[ti]) > 72000:
                n_peaks = count_peaks(rpos_arr)
                if n_peaks > 1:
                    res = fit_gamma(rpos_arr)
                    a1, b1, a2, b2, _, _ = res.x
                    if a1 < a2:
                        gpdf_ti = gpdf(alpha=a1, beta=b1)
                        pos_dist[ti] = gpdf_ti
                    else:
                        gpdf_ti = gpdf(alpha=a2, beta=b2)
                        pos_dist[ti] = gpdf_ti
                    dpeak_fh.write(tname + "\t" + \
                        str(gpdf_ti.alpha) + "\t" + str(gpdf_ti.beta) + "\n")
                    continue
            else:
                dip_stat, pval = diptest(rpos_arr)
                if pval < 0.05 and dip_stat > 0.05:
                    res = fit_gamma(rpos_arr)
                    a1, b1, a2, b2, _, _ = res.x
                    if a1 < a2:
                        gpdf_ti = gpdf(alpha=a1, beta=b1)
                        pos_dist[ti] = gpdf_ti
                    else:
                        gpdf_ti = gpdf(alpha=a2, beta=b2)
                        pos_dist[ti] = gpdf_ti
                    dpeak_fh.write(tname + "\t" + \
                        str(gpdf_ti.alpha) + "\t" + str(gpdf_ti.beta) + "\n")
        dpeak_fh.close()
    else:
        pos_dist = None
        mode_tbl = None
    
    # output the unmapped
    unm_fn = os.path.join(args.out_dir, "unmapped.txt")
    unm_fh = open(unm_fn, 'w')
    for qname in unmapped:
        unm_fh.write(qname + "\n")
    unm_fh.close()

    ti_map_fn = os.path.join(args.out_dir, "ti.pkl")
    with open(ti_map_fn, "wb") as f:
        pickle.dump(ti_map, f)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} building a compatibility matrix")
    cmpt_mat, sur_ctr = build_cmpt_mat(aln_mat, qi_map, pri_lst, tlen_lst, \
                            args.filter, args.five_prime, args.three_prime, \
                            args.surrender, args.target_cover, pos_dist, \
                            args.model)
    
    cmpt_fn = os.path.join(args.out_dir, "cmpt_mat.tsv")
    cmpt_fh = open(cmpt_fn, 'w')
    cmpt_fh.write("#" + str(len(cmpt_mat) - sur_ctr) + "\n")
    score_fn = os.path.join(args.out_dir, "scores.tsv")
    score_fh = open(score_fn, 'w')
    score_fh.write("#" + str(len(cmpt_mat) - sur_ctr) + "\n")
    sur_fn = os.path.join(args.out_dir, "unassigned.txt")
    sur_fh = open(sur_fn, 'w')
    qnames = sorted(qi_map, key=lambda k: qi_map[k])
    tnames = sorted(ti_map, key=lambda k: ti_map[k])
    sur_lines = []
    cmpt_lines = []
    score_lines = []

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out the matrix and scores")
    for qi in tqdm(range(len(qi_map))):
        qname = qnames[qi]
        
        if len(cmpt_mat[qi]) == 0:
            sur_lines.append(qname)
            continue

        cmpt_line = qname
        score_line = qname
        
        for ti in cmpt_mat[qi]:
            cmpt_line += "\t" + tnames[ti]
            score_line += "\t(" + tnames[ti] + "," + str(cmpt_mat[qi][ti]) + ")"
        
        cmpt_lines.append(cmpt_line)
        score_lines.append(score_line)

    sur_fh.write("\n".join(sur_lines) + "\n")
    cmpt_fh.write("\n".join(cmpt_lines) + "\n")
    score_fh.write("\n".join(score_lines) + "\n")
    sur_fh.close()
    cmpt_fh.close()
    score_fh.close()


if __name__ == "__main__":
    main()