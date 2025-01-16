#!/usr/bin/env python

import os
import json
from datetime import datetime
import pickle
import sys
from transigner.commons import RED, GREEN, RESET
from tqdm import tqdm
import scipy.stats as st
import math
import numpy as np
import pandas as pd

def load_cnt(fn):
    cnt_tbl = dict()
    fh = open(fn, 'r')
    for ln in fh:
        clean_ln = ln.strip()
        tmp = clean_ln.split("\t")
        tid = tmp[0]
        cnt = float(tmp[1])
        cnt_tbl[tid] = cnt
    fh.close()
    return cnt_tbl
    
def load_scores(fn, ti_map):
    qi_map = dict()
    ti_set = set()
    first = True
    qi = -1
    with open(fn, 'r') as f:
        for ln in f:
            if first:
                qi_size = int(ln.replace('#', ''))
                cmpt_mat = [dict() for _ in range(qi_size)]
                first = False
                print(datetime.now(), f"{GREEN}PROGRESS{RESET} number of reads:", qi_size)
                continue

            fields = ln.strip().split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qi += 1
                    qname = fields[i]
                    qi_map[qi] = qname # different from prefilter
                else:
                    tname = fields[i].split(',')[0].replace('(', '')
                    score = float(fields[i].split(",")[1].strip().replace(')', ''))
                    ti = ti_map[tname]
                    cmpt_mat[qi][ti] = score
                    ti_set.add(ti)
    # sanity check
    assert qi_size == len(qi_map)
    return cmpt_mat, qi_map, qi_size, ti_set

def load_scores_new(fn, ti_map):
    qi_map = dict()
    ti_set = set()
    first = True
    qi = -1
    prev_qname = None
    with open(fn, 'r') as f:
        for ln in f:
            if first:
                qi_size = int(ln.replace('#', ''))
                cmpt_mat = [dict() for _ in range(qi_size)]
                first = False
                print(datetime.now(), f"{GREEN}PROGRESS{RESET} number of reads:", qi_size)
                continue

            fields = ln.strip().split(",")
            qname = fields[0]
            if prev_qname != qname:
                qi += 1
                qi_map[qi] = qname
                tname = fields[1]
                score = float(fields[2])
                ti = ti_map[tname]
                cmpt_mat[qi][ti] = score
                ti_set.add(ti)
                prev_qname = qname
            else:
                tname = fields[1]
                score = float(fields[2])
                ti = ti_map[tname]
                # print(f'{qname}\t{tname}\t{ti}\t{score}')
                cmpt_mat[qi][ti] = score
                ti_set.add(ti)

    # sanity check
    assert qi_size == len(qi_map)
    return cmpt_mat, qi_map, qi_size, ti_set

def init(qi_size, ti_set):
    rho = dict()
    alpha = [dict() for _ in range(qi_size)]
    ti_size = len(ti_set)
    for ti in ti_set:
        rho[ti] = qi_size / ti_size
    return rho, alpha

def step_e(alpha, rho, qi_size, cmpt_mat, is_naive):
    for qi in range(qi_size):
        alpha_qi = dict()
        loc_ti_set = cmpt_mat[qi]
        for ti in loc_ti_set:
            score = cmpt_mat[qi][ti]
            rho_ti = rho[ti]
            if is_naive:
                alpha_qi[ti] = rho_ti
            else:
                alpha_qi[ti] = score * rho_ti
        z = sum(alpha_qi.values())
        if z == 0:
            for ti in loc_ti_set:
                alpha[qi][ti] = 0
        else:
            for ti in loc_ti_set:
                alpha[qi][ti] = alpha_qi[ti] / z

# output read counts
def step_m(alpha, ti_set, qi_size):
    new_rho = dict()
    for ti in ti_set:
        new_rho[ti] = 0
    for qi in range(qi_size):
        for ti in alpha[qi]:
            new_rho[ti] += alpha[qi][ti]
    return new_rho

def has_converged(old_rho, new_rho, thres):
    delt_rho = 0
    converged = False
    for ti in old_rho:
        # delt_rho += abs(old_rho[ti] - new_rho[ti]
        delt_rho = max(abs(old_rho[ti] - new_rho[ti]), delt_rho)
    if delt_rho < thres:
        converged = True
    return delt_rho, converged
   
def calc_corr(gt_cnt_lst, est_cnt_lst, is_linear):
    if is_linear:
        res = st.pearsonr(gt_cnt_lst, est_cnt_lst)
    else:
        res = st.spearmanr(gt_cnt_lst, est_cnt_lst)
    return res
     
def drop_scores(cmpt_mat, alpha, qi_size, df):
    for qi in range(qi_size):
        n_qi = sum(1 for s in cmpt_mat[qi].values() if s > 0)
        if n_qi <= 1:
            continue
        sigma_qi = 1 / n_qi + (1 / n_qi * df)
        max_alpha = max(cmpt_mat[qi].values())
        max_tis = [k for k, v in cmpt_mat[qi].items() if v == max_alpha]
        ctr = 0
        for ti in cmpt_mat[qi]:
            rf = alpha[qi][ti]
            if rf < sigma_qi:
                cmpt_mat[qi][ti] = 0.0 # no fraction of qi assigned to ti
            else:
                ctr += 1
        if ctr == 0:
            for ti in max_tis:
                cmpt_mat[qi][ti] = max_alpha

def format_rho(rho, tnames):
    est_cnt_tbl = dict()
    for ti in rho:
        tname = tnames[ti].split(".")[0]
        est_cnt_tbl[tname] = rho[ti]
    return est_cnt_tbl

def extract_cnt(gt_cnt_tbl, est_cnt_tbl, take_log):
    gt_cnt_lst = list()
    est_cnt_lst = list()
    eval_tids = list()
    for tid in gt_cnt_tbl:
        gt_cnt = gt_cnt_tbl[tid]
        if gt_cnt < 2: # TODO: is this correct?
            continue
        if tid not in est_cnt_tbl:
            est_cnt = 0
        else:
            est_cnt = est_cnt_tbl[tid]
        eval_tids.append(tid)
        if take_log:
            gt_cnt = np.log2(gt_cnt + 1)
            est_cnt = np.log2(est_cnt + 1)
        gt_cnt_lst.append(gt_cnt)
        est_cnt_lst.append(est_cnt)
    return gt_cnt_lst, est_cnt_lst, eval_tids

def calc_rmse(gt_cnt_lst, est_cnt_lst):
    running_sum = 0
    for gt, est in zip(gt_cnt_lst, est_cnt_lst):
        running_sum += (gt - est) ** 2
    rmse = math.sqrt(running_sum / len(gt_cnt_lst))
    return rmse

def main(args):
    if args.dev:
        if not args.true:
            raise Exception(f"Provide the true read counts in DEV mode")
        gt_cnt_tbl = load_cnt(args.true)

    cmd_fn = os.path.join(args.out_dir, "em_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading target index")
    with open(args.index, 'rb') as f:
        ti_map = pickle.load(f)
    # with open('tmap.csv', 'w') as f:
    #     for tname in ti_map:
    #         f.write(f'{tname},{ti_map[tname]}\n')
    
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading compatibility scores")
    #cmpt_mat, qi_map, qi_size, ti_set = load_scores(args.scores, ti_map)
    cmpt_mat, qi_map, qi_size, ti_set = load_scores_new(args.scores, ti_map)

    tnames = sorted(ti_map, key=lambda k: ti_map[k])

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} starting EM")
    rho, alpha = init(qi_size, ti_set)
    num_iter = 1
    rho_converged = False

    if args.dev:
        log_fn = os.path.join(args.out_dir, 'train.log')
        log_fh = open(log_fn, 'w')
        log_fh.write(f'num_iter,spearman,pearson,rmse\n')

    while num_iter <= args.max_iter:
        step_e(alpha, rho, qi_size, cmpt_mat, args.naive)
        if args.drop and num_iter == 1:
            drop_scores(cmpt_mat, alpha, qi_size, args.drop_fac)
            step_e(alpha, rho, qi_size, cmpt_mat, args.naive)
        new_rho = step_m(alpha, ti_set, qi_size)
        rho = new_rho
        if args.dev:
            est_cnt_tbl = format_rho(rho, tnames)
            gt_cnt_lst, est_cnt_lst, _ = extract_cnt(gt_cnt_tbl, est_cnt_tbl, True) # log
            scc_res = calc_corr(gt_cnt_lst, est_cnt_lst, False)
            gt_cnt_lst, est_cnt_lst, _ = extract_cnt(gt_cnt_tbl, est_cnt_tbl, False)
            pcc_res = calc_corr(gt_cnt_lst, est_cnt_lst, True)
            rmse = calc_rmse(gt_cnt_lst, est_cnt_lst)
            if args.verbose:
                print(f'{num_iter}\t{scc_res.statistic}\t{pcc_res.statistic}\t{rmse}')
            log_fh.write(f'{num_iter},{scc_res.statistic},{pcc_res.statistic},{rmse}\n')

        # delta_rho, converged = has_converged(rho, new_rho, args.rho_thres)
        # if converged:
        #     print(datetime.now(), f"{GREEN}PROGRESS{RESET} converged")
        #     rho_converged = True
        #     break
        num_iter += 1

    if not rho_converged:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} max iter reached")
    
    summary_fn = os.path.join(args.out_dir, "em_summary.txt")
    with open(summary_fn, 'w') as f:
        f.write(f"number of EM iterations: {num_iter - 1}\n")
        # f.write(f"last cumulative rho delta: {delta_rho}\n")
    
    alpha_fn = os.path.join(args.out_dir, "assignments.tsv")
    un_asgn_fn = os.path.join(args.out_dir, "unassigned.txt")
    un_asgn_lines = []
    alpha_lines = []

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out assignments")
    for qi in tqdm(range(qi_size)):
        qname = qi_map[qi]
        if sum(alpha[qi].values()) == 0:
            # print(datetime.now(), f"{RED}WARNING{RESET} read {qname} unassigned")
            un_asgn_lines.append(f"{qname}\n")
            continue
        line = [qname]
        for ti, alpha_ti in alpha[qi].items():
            if alpha_ti > 0:
                tname = tnames[ti]
                line.append(f"\t({tname}, {alpha_ti})")
        alpha_lines.append("\t".join(line) + "\n")
    
    with open(alpha_fn, 'w') as alpha_fh:
        alpha_fh.writelines(alpha_lines)
    
    with open(un_asgn_fn, 'w') as un_asgn_fh:
        un_asgn_fh.writelines(un_asgn_lines)
    
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out abundances")
    rho_fn = os.path.join(args.out_dir, "abundances.csv")
    rho_lines = []
    tnames_set = set(tnames)
    for ti in ti_set:
        tname = tnames[ti]
        rc_ti = rho[ti]
        rho_lines.append(f"{tname},{rc_ti},{rc_ti / qi_size}\n")

    for tname in ti_map:
        if tname not in tnames_set:
            rho_lines.append(f"{tname}\t0.0\t0.0\n")
    
    with open(rho_fn, 'w') as rho_fh:
        rho_fh.write(f'transcript_id,read_count,relative_abundance\n')
        rho_fh.writelines(rho_lines)
        
if __name__ == "__main__":
    main()
