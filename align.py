#!/usr/bin/env python

import argparse
from subprocess import call
from datetime import datetime
from pyGANlib import line, utils
import sys
import os
import json
from commons import RED, GREEN, RESET

def align(args, sN):
    out_fn = os.path.join(args.out_dir, "aln.bam")
    if args.mm2 is not None:
        mm2_cmd = "minimap2 " + args.mm2.strip() + " -N " + str(sN) + " -t " + str(args.threads) \
                    + " " + args.target + " " + args.query + " | samtools sort -@ " + str(args.threads) \
                    + " -o " + out_fn
    else:
        mm2_cmd = "minimap2 -ax map-ont --eqx -N " + str(sN) + " -t " + str(args.threads) \
                    + " " + args.target + " " + args.query + " | samtools sort -@ " + str(args.threads) \
                    + " -o " + out_fn
    print(mm2_cmd)
    call(mm2_cmd, shell=True)
    index_cmd = "samtools index " + out_fn + " -@ " + str(args.threads)
    print(index_cmd)
    call(index_cmd, shell=True)

def calc_max_iso(fn, format, parent_key="gene_id"):
    gene_tbl = dict()
    with open(fn, 'r') as f:
        for ln in f:
            if len(ln.split("\t")) != 9:
                continue
            ln_obj = line.Line(ln, format)
            if ln_obj.feature == "transcript":
                try:
                    assert parent_key in ln_obj.attributes
                except:
                    print(datetime.now(), f"{RED}ERROR{RESET} wrong parent key for transcript features")
                    sys.exit(-1)
                gene_id = ln_obj.attributes[parent_key]
                if gene_id in gene_tbl:
                    gene_tbl[gene_id] += 1
                else:
                    gene_tbl[gene_id] = 1
    max_iso_n = -1
    max_iso_gene = None
    for gene_id in gene_tbl:
        if gene_tbl[gene_id] > max_iso_n:
            max_iso_gene = gene_id
            max_iso_n = gene_tbl[gene_id]
    return max_iso_gene, max_iso_n

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-q', '--query', type=str, help="FASTQ", required=True)
    parser.add_argument('-t', '--target', type=str, help="FASTA", required=True)
    parser.add_argument('-a', '--annot', type=str, help="GTF/GFF", required=False, \
                        default=None)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=True)
    parser.add_argument('-sN', '--sec-num', type=int, help="", \
                        default=181, required=False)
    parser.add_argument('-pad', '--padding', type=int, help="", default=50, \
                        required=False)
    parser.add_argument('-th', '--threads', type=int, help="", default=1, \
                        required=False)
    parser.add_argument('-mm2', type=str, help="", default=None, required=False)
    parser.add_argument('-v', '--verbose', default=False, help="", \
                        required=False, action='store_true')

    args = parser.parse_args()

    cmd_fn = os.path.join(args.out_dir, "align_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    if args.annot is not None:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} calculating max secondary alignments param; -sN will be ignored")
        if args.annot.endswith(("gff", "gff3")):
            format = "gff"
        elif args.annot.endswith("gtf"):
            format = "gtf"
        else:
            print(datetime.now(), f"{RED}ERROR{RESET} unrecognized annotation format")
            sys.exit(-1)
        max_iso_gene, max_iso_n = calc_max_iso(args.annot, format)
        if args.verbose:
            print(f"maximum number of isoforms at {max_iso_gene} locus: {max_iso_n}")
            sN = max_iso_n + args.padding
    else:
        sN = args.sec_num
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} aligning query reads to the target transcriptome")
    # align(args, sN)

if __name__ == "__main__":
    main()
