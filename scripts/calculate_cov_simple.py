#!/usr/bin/env python

# let's keep everything simple

import pysam
import argparse
from tqdm import tqdm
import pyfastx
import sys


def compute_tlen(fn):
    fa = pyfastx.Fasta(fn, build_index=False)
    tlen_tbl = dict()
    tid_set = set()
    for name, seq in fa:
        tlen_tbl[name] = len(seq)
        tid_set.add(name)
    return tlen_tbl, tid_set


def load_aln(fn):
    aln_rlen_tbl = dict()
    aln_rpos_tbl = dict()
    with pysam.AlignmentFile(fn, 'rb') as fh:
        for brec in tqdm(fh):
            rname = brec.query_name
            if brec.is_unmapped:
                continue
            else:
                rst = brec.reference_start
                ren = brec.reference_end
                tname = brec.reference_name
                if rname not in aln_rpos_tbl:
                    # aln_rlen_tbl[rname] = dict()
                    aln_rpos_tbl[rname] = dict()

                # experimental stub
                if tname in aln_rpos_tbl[rname]:
                    aln_rlen_tbl[rname][tname] = max(aln_rlen_tbl[rname][tname], rlen)
                    # aln_rpos_tbl[rname][tname] = (min(aln_rpos_tbl[rname][tname][0], rst),
                    #                               max(aln_rpos_tbl[rname][tname][1], ren))
                else:
                    # aln_rlen_tbl[rname][tname] = rlen
                    aln_rpos_tbl[rname][tname] = (rst, ren)

                if rname == "72956e21-48a8-4315-bfbd-6b4d7a8a67b1":
                    print(aln_rpos_tbl[rname])
    for rname in aln_rpos_tbl:
        aln_rlen_tbl[rname] = dict()
        # print(aln_rpos_tbl[rname])
        for tname in aln_rpos_tbl[rname]:
            rlen = aln_rpos_tbl[rname][tname][1] - aln_rpos_tbl[rname][tname][0]
            aln_rlen_tbl[rname][tname] = rlen
    return aln_rlen_tbl


def load_asgn_top(fn):
    fh = open(fn, 'r')
    asgn_tbl = dict()
    for ln in tqdm(fh):
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        for i in range(len(fields)):
            if i == 0:
                rname = fields[i]
            else:
                tmp = fields[i].split(',')
                tname = tmp[0].strip().replace('(', '')
                ratio = tmp[1].strip().replace(')', '')
                if rname not in asgn_tbl:
                    asgn_tbl[rname] = dict()
                    asgn_tbl[rname][tname] = float(ratio)
                else:
                    asgn_tbl[rname][tname] = float(ratio)
        asgn_tbl[rname] = dict(sorted(asgn_tbl[rname].items(), key=lambda x: x[1], reverse=True))
    fh.close()
    return asgn_tbl


def load_asgn(fn):
    fh = open(fn, 'r')
    asgn_tbl = dict()
    for ln in tqdm(fh):
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        for i in range(len(fields)):
            if i == 0:
                rname = fields[i]
            else:
                tmp = fields[i].split(',')
                tname = tmp[0].strip().replace('(', '')
                ratio = tmp[1].strip().replace(')', '')
                if rname not in asgn_tbl:
                    asgn_tbl[rname] = dict()
                    asgn_tbl[rname][tname] = float(ratio)
                else:
                    asgn_tbl[rname][tname] = float(ratio)
    fh.close()
    return asgn_tbl


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-aln', type=str, help="", required=True)
    parser.add_argument('-asgn', type=str, help="", required=True)
    parser.add_argument('-tgt', type=str, help="", required=True)
    parser.add_argument('-o', '--out-file', type=str, help="", required=False, default="simple_cov.tsv")
    parser.add_argument('-N', '--top-n', type=int, help="", required=False, default=None)
    parser.add_argument('-f', '--frac', type=float, help="", required=False, default=None)
    parser.add_argument('--pad', default=False, help="", required=False, action='store_true')

    args = parser.parse_args()
    print("loading alignments...")
    aln_rlen_tbl = load_aln(args.aln)
    print("loading assignments...")
    if args.top_n is not None:
        print("only considering top N assignments...")
        print("N: " + str(args.top_n))
        asgn_tbl = load_asgn_top(args.asgn)
    elif args.frac is not None:
        print("threshold on alpha: " + str(args.frac))
        asgn_tbl = load_asgn(args.asgn)
    else:
        asgn_tbl = load_asgn(args.asgn)
    print("loading target transcriptome...")
    tlen_tbl, tid_set = compute_tlen(args.tgt)

    bases_tbl = dict()
    for rname in tqdm(asgn_tbl):
        if args.top_n is not None:
            asgn = dict(list(asgn_tbl[rname].items())[:args.top_n])
            norm_fac = sum(asgn.values())
            asgn = {key: value / norm_fac for key, value in asgn.items()}
        elif args.frac is not None:
            asgn = asgn_tbl[rname]
            for k in asgn:
                v = asgn_tbl[rname][k]
                if v >= args.frac:
                    asgn = {k: 1.0}
                    break
        else:
            asgn = asgn_tbl[rname]
        for tname in asgn:
            alpha = asgn_tbl[rname][tname]
            try:
                aln_rlen = aln_rlen_tbl[rname][tname]
            except:
                print(aln_rlen_tbl[rname])
                print(rname)
                print(tname)
            num_b = alpha * aln_rlen
            if tname not in bases_tbl:
                bases_tbl[tname] = num_b
            else:
                bases_tbl[tname] += num_b

    cov_tbl = dict()
    for tname in bases_tbl:
        tlen = tlen_tbl[tname]
        cov_tbl[tname] = bases_tbl[tname] / tlen

    out_fh = open(args.out_file, 'w')

    tid_curr_set = set()
    for tname in cov_tbl:
        tid_curr_set.add(tname)
        out_fh.write(tname + "\t" + str(cov_tbl[tname]) + "\n")
    if args.pad:
        remain_tid_set = tid_set - tid_curr_set
        for tid in remain_tid_set:
            out_fh.write(tid + "\t0\n")

    out_fh.close()


if __name__ == "__main__":
    main()
