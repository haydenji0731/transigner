#!/usr/bin/env python

import sys


def load_gt(fn):
    fh = open(fn, 'r')
    gt_cnt_tbl = dict()
    for line in fh:
        clean_ln = line.strip()
        tid = clean_ln.split("\t")[0]
        cnt = int(clean_ln.split("\t")[1])
        gt_cnt_tbl[tid] = cnt
    fh.close()
    return gt_cnt_tbl


def load_qry(fn):
    fh = open(fn, 'r')
    qry_cnt_tbl = dict()
    for line in fh:
        clean_ln = line.strip()
        tid = clean_ln.split("\t")[0]
        cnt = float(clean_ln.split("\t")[1])
        qry_cnt_tbl[tid] = cnt
    fh.close()
    return qry_cnt_tbl


def main(gt_fn, qry_fn, out_fn):
    gt_cnt_tbl = load_gt(gt_fn)
    qry_cnt_tbl = load_qry(qry_fn)
    out_fh = open(out_fn, 'w')
    qry_set = set(qry_cnt_tbl)
    gt_set = set(gt_cnt_tbl)
    for tid in gt_cnt_tbl:
        if tid not in qry_cnt_tbl:
            out_fh.write(tid + "\t" + str(0) + "\n")
        else:
            cnt = qry_cnt_tbl[tid]
            out_fh.write(tid + "\t" + str(cnt) + "\n")
    diff_set = qry_set.difference(gt_set)
    for tid in diff_set:
        cnt = qry_cnt_tbl[tid]
        out_fh.write(tid + "\t" + str(cnt) + "\n")
    out_fh.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
