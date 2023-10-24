#!/usr/bin/env python

import sys


def load_gt(fn):
    gt_cnt_tbl = dict()
    fh = open(fn, 'r')
    tot = 0
    for line in fh:
        clean_ln = line.strip()
        tid = clean_ln.split("\t")[0]
        cnt = int(clean_ln.split("\t")[1])
        if cnt == 0:
            print("ERROR: zero count ground truth tx is not expected")
            sys.exit()
        gt_cnt_tbl[tid] = cnt
        tot += cnt
    return gt_cnt_tbl


def load_qry(fn):
    qry_cnt_tbl = dict()
    fh = open(fn, 'r')
    tot_cnt = 0
    for line in fh:
        clean_ln = line.strip()
        tid = clean_ln.split("\t")[0]
        cnt = float(clean_ln.split("\t")[1])
        qry_cnt_tbl[tid] = cnt
        tot_cnt += cnt
    print("total query count: " + str(tot_cnt))
    return qry_cnt_tbl


def compute_loss(gt_cnt_tbl, qry_cnt_tbl, gt_id_set=None):
    tae = 0
    ctr = 0
    tot = 0
    ignored = 0
    ignored_cnt = 0
    for tid in qry_cnt_tbl:
        if gt_id_set:
            if tid not in gt_id_set:
                # ignored += 1
                cnt = qry_cnt_tbl[tid]
                if cnt > 0:
                    ignored += 1
                ignored_cnt += cnt
                continue
        ctr += 1
        if tid in gt_cnt_tbl:
            tot += gt_cnt_tbl[tid]
            e = abs(gt_cnt_tbl[tid] - qry_cnt_tbl[tid])
        else:
            e = qry_cnt_tbl[tid]
        tae += e
    mae = tae / len(gt_cnt_tbl)
    print("number of ignored transcripts: " + str(ignored))
    print("total count assigned to these ignored transcripts: " + str(ignored_cnt))
    return tae, mae


def load_id(fn):
    fh = open(fn, 'r')
    id_set = set()
    for line in fh:
        clean_ln = line.strip()
        id_set.add(clean_ln)
    fh.close()
    return id_set


def main(gt_fn, gt_id_fn, qry_fn):
    gt_cnt_tbl = load_gt(gt_fn)
    qry_cnt_tbl = load_qry(qry_fn)
    gt_id_set = load_id(gt_id_fn)
    # sanity check
    for tid in gt_cnt_tbl:
        if tid not in qry_cnt_tbl:
            print("FATAL ERROR: all ground truth txs must be present in the query count table...")
            sys.exit()
    tae, mae = compute_loss(gt_cnt_tbl, qry_cnt_tbl, gt_id_set)
    print(str(tae) + "\t" + str(mae))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
