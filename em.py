#!/usr/bin/env python

from prefilter import load_cmpt_tbl, load_score_tbl
import pickle
import argparse
from tqdm import tqdm
import os


def step_e(assignment, abundance, score_tbl, cmpt_tbl, reads_lst, use_score):
    reads_num = len(reads_lst)
    for i in tqdm(range(reads_num)):
        tmp = dict()
        ti_set = cmpt_tbl[i]
        for ti in ti_set:
            score = score_tbl[i][ti]
            rho = abundance[ti]
            if use_score:
                tmp[ti] = score * rho
            else:
                tmp[ti] = rho
        z = sum(tmp.values())
        if z == 0:
            print(reads_lst[i])
        for ti in ti_set:
            assignment[i][ti] = tmp[ti] / z


def step_m(assignment, cmpt_tbl, reads_lst, glob_ti_set):
    updated_abundance = dict()
    for ti in glob_ti_set:
        updated_abundance[ti] = 0

    reads_num = len(reads_lst)
    for i in tqdm(range(reads_num)):
        ti_set = cmpt_tbl[i]
        for ti in ti_set:
            alpha = assignment[i][ti]
            updated_abundance[ti] += alpha

    for ti in glob_ti_set:
        # reads_num = total # of aligned reads
        updated_abundance[ti] /= reads_num
    return updated_abundance


def has_converged(old_abundance, new_abundance, thres, glob_ti_set):
    delt = 0
    converged = False
    for ti in glob_ti_set:
        delt += abs(old_abundance[ti] - new_abundance[ti])
    if delt < thres:
        converged = True
    return delt, converged


def drop_scores(score_tbl, assignment, reads_lst):
    reads_num = len(reads_lst)
    for i in tqdm(range(reads_num)):
        score_d = score_tbl[i]
        # TODO: did this fix?
        cmpt_num = sum(1 for j in score_d.values() if j != 0)
        thres = 1 / cmpt_num
        for ti in score_d:
            read_frac = assignment[i][ti] + 0.0000000000000001
            if read_frac < thres:
                score_tbl[i][ti] = 0


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-cmpt', '--cmpt-tsv', type=str, help="", required=True)
    parser.add_argument('-scores', '--scores-tsv', type=str, help="", required=True)
    parser.add_argument('-i', '--index', type=str, help="", required=True)
    parser.add_argument('-thres', type=float, help="minimum cumulative rho change", required=False, default=0.005)
    parser.add_argument('-max-iter', type=int, help="maximum number of EM iterations", required=False, default=100)
    parser.add_argument('-o', '--out_dir', type=str, help="", required=True)
    parser.add_argument('--use-score', default=False, help="", required=False, action='store_true')
    parser.add_argument('--drop', default=False, help="", required=False, action='store_true')
    parser.add_argument('-l', '--len', type=int, help="", required=True)

    args = parser.parse_args()

    print("loading transcriptome index")
    with open(args.index, "rb") as f:
        tx_tbl = pickle.load(f)
    tmp_lst = list()
    for tname in tx_tbl:
        tmp_lst.append((tname, tx_tbl[tname]))
    sorted_tmp_lst = sorted(tmp_lst, key=lambda x: x[1])
    tx_lst = [x for (x, i) in sorted_tmp_lst]

    print("loading compatibility matrix")
    cmpt_tbl, reads_lst = load_cmpt_tbl(args.cmpt_tsv, args.len, tx_tbl)

    print("loading score matrix")
    score_tbl, reads_lst_2 = load_score_tbl(args.scores_tsv, args.len, tx_tbl)

    # sanity check
    assert len(reads_lst) == len(reads_lst_2)

    print("initializing...")
    abundance = dict()
    assignment = [dict() for _ in reads_lst]
    reads_num = len(reads_lst)
    glob_ti_set = set()
    for i in tqdm(range(reads_num)):
        ti_set = cmpt_tbl[i]
        for ti in ti_set:
            glob_ti_set.add(ti)
    glob_ti_num = len(glob_ti_set)
    for ti in glob_ti_set:
        abundance[ti] = 1 / glob_ti_num

    num_iter = 1
    print("starting EM...")
    rho_conv = False
    while num_iter <= args.max_iter:
        print("Iteration #: " + str(num_iter))
        print("E step\n-------------------------------")
        step_e(assignment, abundance, score_tbl, cmpt_tbl, reads_lst, args.use_score)
        if args.drop and num_iter == 1:
            print("Dropping some low-compatibility Txs...")
            drop_scores(score_tbl, assignment, reads_lst)
            step_e(assignment, abundance, score_tbl, cmpt_tbl, reads_lst, args.use_score)
        print("M step\n-------------------------------")
        new_abundance = step_m(assignment, cmpt_tbl, reads_lst, glob_ti_set)
        delta, converged = has_converged(abundance, new_abundance, args.thres, glob_ti_set)
        print("Rho Delta: " + str(delta))
        abundance = new_abundance
        if converged:
            print("Converged w.r.t. minimum cumulative rho delta criteria")
            rho_conv = True
            break
        num_iter += 1
    if not rho_conv:
        print("Converged w.r.t. maximum # of iterations criteria")

    print("writing out relative abundances")
    abundance_fn = os.path.join(args.out_dir, "abundances.tsv")
    abundance_fh = open(abundance_fn, 'w')
    for ti in tqdm(glob_ti_set):
        tname = tx_lst[ti]
        rho = abundance[ti]
        cnt = rho * reads_num
        abundance_fh.write(tname + "\t" + str(rho) + "\t" + str(cnt) + "\n")
    abundance_fh.close()

    print("writing out read-to-tx assignments")
    assignment_fn = os.path.join(args.out_dir, "assignments.tsv")
    assignment_fh = open(assignment_fn, 'w')
    max_fn = os.path.join(args.out_dir, "max_assignments.tsv")
    max_fh = open(max_fn, 'w')
    for i in tqdm(range(reads_num)):
        qname = reads_lst[i]
        ti_set = cmpt_tbl[i]
        assignment_fh.write(qname)
        max_fh.write(qname)
        tmp_lst = list()
        for ti in ti_set:
            tname = tx_lst[ti]
            alpha = assignment[i][ti]
            tmp_lst.append((tname, alpha))
            if alpha > 0:
                assignment_fh.write("\t(" + tname + ", " + str(alpha) + ")")
        sorted_tmp_lst = sorted(tmp_lst, key=lambda x: x[1], reverse=True)
        max_tname, max_alpha = sorted_tmp_lst[0]
        max_fh.write("\t(" + max_tname + ", " + str(max_alpha) + ")\n")
        assignment_fh.write("\n")
    assignment_fh.close()


if __name__ == "__main__":
    main()













