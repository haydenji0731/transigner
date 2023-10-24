#!/usr/bin/env python

import argparse
import sys
from format import Tool
from enum import Enum
from tqdm import tqdm
import re


class Mode(Enum):
    # no penalty to unassembled
    MODE_A = 1
    # yes penalty to unassembled
    MODE_B = 2

    def __str__(self):
        return self.name


def load_reads(fn):
    read_set = set()
    fh = open(fn, 'r')
    for line in tqdm(fh):
        read = line.strip()
        read_set.add(read)
    # print("total # of reads: " + str(len(read_set)))
    return read_set


def load_full_gtf(fn):
    full_rid_set = set()
    fh = open(fn, 'r')
    for line in fh:
        if line[0] == '#':
            continue
        clean_ln = line.strip()
        fields = clean_ln.split("\t")
        feature = fields[2]
        infos = fields[8].split(";")
        if feature == "transcript":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    rid = val.split(".")[0]
                    full_rid_set.add(rid)
                    break
    fh.close()
    return full_rid_set


def load_asgn_full(asgn_fn, tool, is_soft):
    asgn_tbl = dict()
    asgn_fh = open(asgn_fn, 'r')
    for line in asgn_fh:
        # TRANSIGNER
        if tool.value == 1:
            fields = line.strip().split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qid = fields[i]
                else:
                    tid = fields[i].split(",")[0].replace('(', '')
                    tid = tid.split(".")[0]
                    score = float(fields[i].split(",")[1].replace(')', ''))
                    if is_soft:
                        if qid not in asgn_tbl:
                            asgn_tbl[qid] = [(tid, score)]
                        else:
                            asgn_tbl[qid].append((tid, score))
                    else:
                        asgn_tbl[qid] = tid
        # BAMBU
        elif tool.value == 3:
            fields = line.strip().split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qid = fields[i]
                else:
                    tid = fields[i]
                    if tid == "None":
                        # asgn_tbl[qid] = None
                        continue
                    tid = tid.split(".")[0]
                    if qid not in asgn_tbl:
                        asgn_tbl[qid] = [tid]
                    else:
                        asgn_tbl[qid].append(tid)
        else:
            print("invalid tool option argument provided")
            sys.exit()
    return asgn_tbl


def load_asgn(fn, tool, is_soft):
    asgn_tbl = dict()
    fh = open(fn, 'r')
    for line in fh:
        # TRANSIGNER
        if tool.value == 1:
            fields = line.strip().split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qid = fields[i]
                else:
                    tid = fields[i].split(",")[0].replace('(', '')
                    score = float(fields[i].split(",")[1].replace(')', ''))
                    if is_soft:
                        if qid not in asgn_tbl:
                            asgn_tbl[qid] = [(tid, score)]
                        else:
                            asgn_tbl[qid].append((tid, score))
                    else:
                        asgn_tbl[qid] = tid
        # BAMBU
        elif tool.value == 3:
            fields = line.strip().split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qid = fields[i]
                else:
                    tid = fields[i]
                    if tid == "None":
                        asgn_tbl[qid] = None
                    else:
                        if qid not in asgn_tbl:
                            asgn_tbl[qid] = [tid]
                        else:
                            asgn_tbl[qid].append(tid)
        # FLAIR
        elif tool.value == 4:
            tid = line.split("\t")[0]
            qid_ln = line.split("\t")[1].strip()
            qid_lst = qid_ln.split(",")
            for qid in qid_lst:
                asgn_tbl[qid] = tid
        else:
            print("invalid tool option argument provided")
            sys.exit()
    return asgn_tbl


def load_ref(fn, is_loose):
    ref_tbl = dict()
    fh = open(fn, 'r')
    first = True
    ref_eq_tbl = dict()
    glob_rid_set = set()
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        rid = line.split("\t")[1]
        if rid not in ref_eq_tbl:
            ref_eq_tbl[rid] = False
        code = line.split("\t")[2].strip()
        if rid == '.':
            ref_tbl[tid] = None
            continue
        if is_loose:
            glob_rid_set.add(rid)
            if code == '=':
                ref_eq_tbl[rid] = True
            ref_tbl[tid] = (rid, code)
        else:
            if code == '=':
                glob_rid_set.add(rid)
                ref_tbl[tid] = rid
            else:
                ref_tbl[tid] = None
    fh.close()
    return ref_tbl, glob_rid_set, ref_eq_tbl


def get_gt_from_qid(qid):
    qid_pieces = qid.split("_")
    if qid_pieces[0] == "unassigned":
        gt = qid_pieces[0] + "_" + qid_pieces[1] + "_" + qid_pieces[2]
    else:
        gt = qid_pieces[0] + "_" + qid_pieces[1]
    return gt


def evaluate_full(asgn_tbl, read_set, full_rid_set, tool, mode, is_soft):
    tp = 0
    fp = 0
    fn = 0
    asgn_qid_lst = list(asgn_tbl.keys())
    for i in tqdm(range(len(asgn_qid_lst))):
        qid = asgn_qid_lst[i]
        gt = get_gt_from_qid(qid)
        # TRANSIGNER
        if tool.value == 1:
            if is_soft:
                tmp = asgn_tbl[qid]
                found = False
                for tid, score in tmp:
                    if gt == tid:
                        tp += score
                        fn += 1 - score
                        found = True
                    else:
                        fp += score
                if not found:
                    if mode.value == 1:
                        if gt in full_rid_set:
                            fn += 1
                    elif mode.value == 2:
                        fn += 1
                    else:
                        print("invalid mode option argument provided")
                        sys.exit()
            else:
                tid = asgn_tbl[qid]
                if tid == gt:
                    tp += 1
                else:
                    fp += 1
                    if mode.value == 1:
                        if gt in full_rid_set:
                            fn += 1
                    elif mode.value == 2:
                        fn += 1
                    else:
                        print("invalid mode option argument provided")
                        sys.exit()
        # BAMBU
        elif tool.value == 3:
            tid_lst = asgn_tbl[qid]
            tid_lst_len = len(tid_lst)
            if gt in tid_lst:
                tp += 1 / tid_lst_len
                fp += (tid_lst_len - 1) / tid_lst_len
                fn += (tid_lst_len - 1) / tid_lst_len
            else:
                fp += 1
                if mode.value == 1:
                    if gt in full_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
                else:
                    print("invalid mode option argument provided")
                    sys.exit()
        else:
            print("invalid tool option argument provided")
            sys.exit()
    un_asgn_ed = read_set.difference(set(asgn_qid_lst))
    if mode.value == 1:
        for qid in un_asgn_ed:
            gt = get_gt_from_qid(qid)
            if gt in full_rid_set:
                fn += 1
    elif mode.value == 2:
        fn += len(un_asgn_ed)
    else:
        print("invalid mode option argument provided")
        sys.exit()
    return tp, fp, fn


def evaluate_strict(asgn_tbl, ref_tbl, glob_rid_set, read_set, tool, mode, is_soft):
    tp = 0
    fp = 0
    fn = 0
    asgn_qid_lst = list(asgn_tbl.keys())
    for i in tqdm(range(len(asgn_qid_lst))):
        qid = asgn_qid_lst[i]
        gt = get_gt_from_qid(qid)
        # TRANSIGNER
        if tool.value == 1:
            if is_soft:
                asgn_lst = asgn_tbl[qid]
                found = False
                tp_portion = 0
                for tid, score in asgn_lst:
                    rid = ref_tbl[tid]
                    if gt == rid:
                        tp_portion += score
                        found = True
                    else:
                        fp += score
                if not found:
                    if mode.value == 1:
                        if gt in glob_rid_set:
                            fn += 1
                    elif mode.value == 2:
                        fn += 1
                    else:
                        print("invalid mode option argument provided")
                        sys.exit()
                else:
                    # TODO: discuss this w/ ela
                    tp += tp_portion
                    fn += (1 - tp_portion)
            # hard, 1-to-1 assignment
            else:
                tid = asgn_tbl[qid]
                rid = ref_tbl[tid]
                if gt == rid:
                    tp += 1
                else:
                    fp += 1
                    if mode.value == 1:
                        if gt in glob_rid_set:
                            fn += 1
                    elif mode.value == 2:
                        fn += 1
                    else:
                        print("invalid mode option argument provided")
                        sys.exit()
        # BAMBU - only soft available
        elif tool.value == 3:
            tid_lst = asgn_tbl[qid]
            if tid_lst is None:
                if mode.value == 1:
                    if gt in glob_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
                else:
                    print("invalid mode option argument provided")
                    sys.exit()
                continue
            tid_lst_len = len(tid_lst)
            found = False
            tp_portion = 0
            for tid in tid_lst:
                rid = ref_tbl[tid]
                if gt == rid:
                    tp_portion += 1 / tid_lst_len
                    found = True
                else:
                    fp += 1 / tid_lst_len
            if not found:
                if mode.value == 1:
                    if gt in glob_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
                else:
                    print("invalid mode option argument provided")
                    sys.exit()
            else:
                tp += tp_portion
                fn += (1 - tp_portion)
        # FLAIR - only hard available
        elif tool.value == 4:
            tid = asgn_tbl[qid]

            # processing FLAIR transcript ids
            tmp = tid.split("_")
            type_1 = False
            type_2 = False
            chr_idx = None

            for i in range(len(tmp)):
                el = tmp[i]
                if el[0:3] == "chr":
                    type_2 = True
                    chr_idx = i
            if not type_2:
                type_1 = True
            if type_1:
                tmp = tmp[:-1]
            else:
                del_idx = -1 * (len(tmp) - chr_idx)
                tmp = tmp[:del_idx]
            tid_new = ""
            for i in range(len(tmp)):
                if i != len(tmp) - 1:
                    tid_new += tmp[i]
                    tid_new += "_"
                else:
                    tid_new += tmp[i]
            tid = tid_new

            rid = ref_tbl[tid]
            if gt == rid:
                tp += 1
            else:
                fp += 1
                if mode.value == 1:
                    if gt in glob_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
    un_asgn_ed = read_set.difference(set(asgn_qid_lst))
    if mode.value == 1:
        for qid in un_asgn_ed:
            gt = get_gt_from_qid(qid)
            if gt in glob_rid_set:
                fn += 1
    elif mode.value == 2:
        fn += len(un_asgn_ed)
    else:
        print("invalid mode option argument provided")
        sys.exit()
    return tp, fp, fn


# TODO: fix this
def evaluate_loose(asgn_tbl, ref_tbl, glob_rid_set, read_set, tool, mode, is_soft, ref_eq_tbl):
    tp = 0
    fp = 0
    fn = 0
    asgn_qid_lst = list(asgn_tbl.keys())
    for i in tqdm(range(len(asgn_qid_lst))):
        qid = asgn_qid_lst[i]
        gt = get_gt_from_qid(qid)

        eq_exists = False
        if gt in ref_eq_tbl:
            if ref_eq_tbl[gt]:
                eq_exists = True

        if tool.value == 1:
            if is_soft:
                tmp = asgn_tbl[qid]
                found = False
                tp_portion = 0
                for tid, score in tmp:
                    if tid in ref_tbl:
                        rid, code = ref_tbl[tid]
                        if gt == rid:
                            if eq_exists:
                                if code == '=':
                                    tp_portion += score
                                    found = True
                                else:
                                    fp += score
                            else:
                                tp_portion += score
                                found = True
                        else:
                            fp += score
                    else:
                        fp += score
                if not found:
                    if mode.value == 1:
                        if gt in glob_rid_set:
                            fn += 1
                    elif mode.value == 2:
                        fn += 1
                else:
                    tp += tp_portion
                    fn += (1 - tp_portion)
        elif tool.value == 3:
            tid_lst = asgn_tbl[qid]
            tid_lst_len = len(tid_lst)
            found = False
            tp_portion = 0
            for tid in tid_lst:
                if tid in ref_tbl:
                    rid, code = ref_tbl[tid]
                    if gt == rid:
                        if eq_exists:
                            if code == '=':
                                tp_portion += 1 / tid_lst_len
                                found = True
                            else:
                                fp += 1 / tid_lst_len
                        else:
                            tp_portion += 1 / tid_lst_len
                            found = True
                    else:
                        fp += 1 / tid_lst_len
                else:
                    fp += 1 / tid_lst_len
            if not found:
                if mode.value == 1:
                    if gt in glob_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
            else:
                tp += tp_portion
                fn += (1 - tp_portion)
        elif tool.value == 4:
            tid = asgn_tbl[qid]
            found = False
            if tid in ref_tbl:
                rid = ref_tbl[tid]
                if gt == rid:
                    if eq_exists:
                        if code == '=':
                            tp += 1
                            found = True
                        else:
                            fp += 1
                else:
                    fp += 1
            else:
                fp += 1
            if not found:
                if mode.value == 1:
                    if gt in glob_rid_set:
                        fn += 1
                elif mode.value == 2:
                    fn += 1
    un_asgn_ed = read_set.difference(set(asgn_qid_lst))
    if mode.value == 1:
        for qid in un_asgn_ed:
            gt = get_gt_from_qid(qid)
            if gt in glob_rid_set:
                fn += 1
    elif mode.value == 2:
        fn += len(un_asgn_ed)
    else:
        print("invalid mode option argument provided")
        sys.exit()
    return tp, fp, fn


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-ref', '--ref-tsv', type=str, help="", required=False, default=None)
    parser.add_argument('-asgn', '--asgn-tsv', type=str, help="", required=True)
    parser.add_argument('--tool', type=lambda tool: Tool[tool], help="", required=True, choices=list(Tool))
    parser.add_argument('--mode', type=lambda mode: Mode[mode], help="", required=True, choices=list(Mode),
                        default=None)
    parser.add_argument('--loose', default=False, help="", required=False, action='store_true')
    parser.add_argument('--full', default=False, help="", required=False, action='store_true')
    parser.add_argument('--full-gtf', default=None, help="", required=False, type=str)
    parser.add_argument('--soft', default=False, help="", required=False, action='store_true')
    parser.add_argument('-rd', '--reads', type=str, help="", required=True)
    parser.add_argument('--silent', default=False, help="", required=False, action='store_true')
    args = parser.parse_args()

    # some argument checks
    # NanoCount doesn't output read-to-tx mappings
    if args.tool.name == "NANOCOUNT":
        print("cannot run an evaluation on NanoCount output")
        sys.exit()

    # full mode
    if args.full:
        if args.full_gtf is None:
            print("ARGUMENT ERROR: please provide a full gtf")
            sys.exit()
        if args.ref_tsv is not None:
            print("de novo-to-reference tx mappings will be ignored when the full annotation was provided for "
                  "read assignment")
        if args.loose:
            print("loose option can't be used for when the full annotation was provided for read assignment")
            sys.exit()

    # tool-specific modes
    if args.tool.name not in ['TRANSIGNER', 'BAMBU']:
        # FLAIR only outputs hard, 1-to-1 mappings
        if args.soft:
            print("soft option is only available for TranSigner or Bambu outputs")
            sys.exit()
        if args.full:
            print("full option is only available for TranSigner or Bambu outputs")
            sys.exit()

    if args.tool.name == 'BAMBU' and not args.soft:
        print("Bambu outputs only soft assignments.")
        sys.exit()

    # load reads
    if not args.silent:
        print("loading query read list...")
    read_set = load_reads(args.reads)

    if args.full:
        if not args.silent:
            print("loading assignments")
        asgn_tbl = load_asgn_full(args.asgn_tsv, args.tool, args.soft)
        if not args.silent:
            print("loading target transcriptome gtf")
        full_rid_set = load_full_gtf(args.full_gtf)
        if not args.silent:
            print("evaluating")
        tp, fp, fn = evaluate_full(asgn_tbl, read_set, full_rid_set, args.tool, args.mode, args.soft)
    else:
        if not args.silent:
            print("loading tx-to-reference mappings")
        ref_tbl, glob_rid_set, ref_eq_tbl = load_ref(args.ref_tsv, args.loose)
        if not args.silent:
            print("loading assignments")
        asgn_tbl = load_asgn(args.asgn_tsv, args.tool, args.soft)
        if not args.silent:
            print("evaluating")
        if args.loose:
            tp, fp, fn = evaluate_loose(asgn_tbl, ref_tbl, glob_rid_set, read_set, args.tool, args.mode, args.soft,
                                        ref_eq_tbl)
        else:
            tp, fp, fn = evaluate_strict(asgn_tbl, ref_tbl, glob_rid_set, read_set, args.tool, args.mode, args.soft)
    precision = tp / (tp + fp)
    sensitivity = tp / (tp + fn)
    f1 = (2 * precision * sensitivity) / (precision + sensitivity)
    print(str(tp) + "\t" + str(fp) + "\t" + str(fn) + "\t" +
          str(sensitivity) + "\t" + str(precision) + "\t" + str(f1))


if __name__ == "__main__":
    main()
