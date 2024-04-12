#!/usr/bin/env python

import sys
import argparse
from enum import Enum


class Tool(Enum):
    TRANSIGNER = 1
    NANOCOUNT = 2
    BAMBU = 3
    FLAIR = 4
    PUSH = 5
    ESPRESSO = 6

    def __str__(self):
        return self.name


def load_ref(fn, is_loose):
    ref_tbl = dict()
    fh = open(fn, 'r')
    first = True
    ref_eq_tbl = dict()
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        rid = line.split("\t")[1]
        if rid not in ref_eq_tbl:
            ref_eq_tbl[rid] = False
        code = line.split("\t")[2].strip()
        # TODO: fix this
        if rid == '.':
            ref_tbl[tid] = None
            continue
        if is_loose:
            if code == '=':
                ref_eq_tbl[rid] = True
            ref_tbl[tid] = (rid, code)
        else:
            if code == '=':
                ref_tbl[tid] = rid
            else:
                ref_tbl[tid] = None
    fh.close()
    return ref_tbl, ref_eq_tbl


def load_cnt(fn, tool, is_full):
    cnt_tbl = dict()
    fh = open(fn, 'r')
    first = True
    tot = 0
    for line in fh:
        tid = line.split("\t")[0]
        if is_full:
            tid = tid.split(".")[0]
        # TRANSIGNER
        if tool.value == 1:
            cnt = float(line.split("\t")[2])
        # POST-PUSH TRANSIGNER
        elif tool.value == 5:
            cnt = float(line.split("\t")[1])
        # NANOCOUNT or BAMBU
        elif tool.value == 2 or tool.value == 3:
            if first:
                first = False
                continue
            cnt = float(line.split("\t")[2])
        # FLAIR
        elif tool.value == 4:
            if first:
                first = False
                continue
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

            cnt = float(line.split("\t")[1])
        elif tool.value == 6:
            if first:
                first = False
                continue
            clean_ln = line.strip()
            cnt = float(clean_ln.split("\t")[3])
        else:
            print("invalid tool option argument provided")
            sys.exit()
        cnt_tbl[tid] = cnt
        tot += cnt
    print("total abundance: " + str(tot))
    return cnt_tbl


def map_cnt2ref(cnt_tbl, ref_tbl):
    mapped_cnt_tbl = dict()
    tot = 0
    for tid in cnt_tbl:
        cnt = cnt_tbl[tid]
        tot += cnt
        rid = ref_tbl[tid]
        if rid is None:
            mapped_cnt_tbl[tid] = cnt
        else:
            if rid in mapped_cnt_tbl:
                # print("duplicate mapping detected")
                # print(rid)
                mapped_cnt_tbl[rid] += cnt
            else:
                mapped_cnt_tbl[rid] = cnt
    return mapped_cnt_tbl


def map_cnt2ref_loose(cnt_tbl, ref_tbl, ref_eq_tbl):
    mapped_cnt_tbl = dict()
    for tid in cnt_tbl:
        tmp = ref_tbl[tid]
        if tmp is None:
            mapped_cnt_tbl[tid] = cnt_tbl[tid]
            continue
        rid, code = tmp
        if ref_eq_tbl[rid]:
            if code == '=':
                if rid not in mapped_cnt_tbl:
                    mapped_cnt_tbl[rid] = cnt_tbl[tid]
                else:
                    mapped_cnt_tbl[rid] += cnt_tbl[tid]
            else:
                mapped_cnt_tbl[tid] = cnt_tbl[tid]
        else:
            if rid not in mapped_cnt_tbl:
                mapped_cnt_tbl[rid] = cnt_tbl[tid]
            else:
                mapped_cnt_tbl[rid] += cnt_tbl[tid]
    return mapped_cnt_tbl


def write_cnt(cnt_tbl, tool, out_fn):
    out_fh = open(out_fn, 'w')
    for tid in cnt_tbl:
        # TRANSIGNER or NANOCOUNT
        if tool.value == 1 or tool.value == 2:
            cnt = cnt_tbl[tid]
            out_fh.write(tid + "\t" + str(cnt) + "\n")
        # BAMBU
        # TODO: fix this (i.e., remove)
        elif tool.value == 3:
            cnt = cnt_tbl[tid]
            out_fh.write(tid + "\t" + str(cnt) + "\n")
        else:
            print("invalid tool option argument provided")
            sys.exit()
    out_fh.close()


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-ref', '--ref-tsv', type=str, help="", required=False, default=None)
    parser.add_argument('-cnt', '--cnt-tsv', type=str, help="", required=True)
    parser.add_argument('--tool', type=lambda tool: Tool[tool], help="", required=True, choices=list(Tool))
    parser.add_argument('--full', default=False, help="", required=False, action='store_true')
    parser.add_argument('--loose', default=False, help="", required=False, action='store_true')
    parser.add_argument('-o', '--output-fn', type=str, help="", required=True)
    args = parser.parse_args()

    if args.tool.name == "FLAIR" and args.full:
        print("full option is not compatible with FLAIR outputs")
        sys.exit()

    if not args.full:
        if args.ref_tsv is None:
            print("formatting count tsv files for incomplete annotation requires a de novo-to-reference tx mappings")
            sys.exit()
        ref_tbl, ref_eq_tbl = load_ref(args.ref_tsv, args.loose)

    if args.full and args.loose:
        print("loose option is not compatible with outputs obtained from full, complete annotation")
        sys.exit()

    cnt_tbl = load_cnt(args.cnt_tsv, args.tool, args.full)

    if not args.full:
        if args.loose:
            mapped_cnt_tbl = map_cnt2ref_loose(cnt_tbl, ref_tbl, ref_eq_tbl)
        else:
            mapped_cnt_tbl = map_cnt2ref(cnt_tbl, ref_tbl)
        out_fh = open(args.output_fn, 'w')
        for tid in mappe:d_cnt_tbl:
            cnt = mapped_cnt_tbl[tid]
            out_fh.write(tid + "\t" + str(cnt) + "\n")
        out_fh.close()
    else:
        write_cnt(cnt_tbl, args.tool, args.output_fn)


if __name__ == "__main__":
    main()
