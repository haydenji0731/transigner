#!/usr/bin/env python

from transigner import align, em, prefilter
import argparse


def parse():
    parser = argparse.ArgumentParser(description="")
    subparsers = parser.add_subparsers(dest='module', \
                                       help="specify one of three transigner modules")
    parser_align = subparsers.add_parser('align', help="arguments for the align module")
    parser_align.add_argument('-q', '--query', type=str, help="FASTQ", required=True)
    parser_align.add_argument('-t', '--target', type=str, help="FASTA", required=True)
    parser_align.add_argument('-a', '--annot', type=str, help="GTF/GFF", required=False, \
                        default=None)
    parser_align.add_argument('-d', '--out-dir', type=str, help="", required=False, default=".")
    parser_align.add_argument('-o', '--out-file', type=str, help="", required=True)
    parser_align.add_argument('-sN', '--sec-num', type=int, help="", \
                        default=181, required=False)
    parser_align.add_argument('-pad', '--padding', type=int, help="", default=50, \
                        required=False)
    parser_align.add_argument('-p', '--threads', type=int, help="", default=1, \
                        required=False)
    parser_align.add_argument('-mm2', type=str, help="", default=None, required=False)
    parser_align.add_argument('-v', '--verbose', default=False, help="", \
                        required=False, action='store_true')
    
    parser_pref = subparsers.add_parser('prefilter', help="arguments for the prefilter module")
    parser_pref.add_argument('-a', '--aln', type=str, help="BAM or SAM", required=True)
    parser_pref.add_argument('-t', '--target', type=str, help="FASTA", required=True)
    parser_pref.add_argument('-o', '--out-dir', type=str, help="", required=True)
    parser_pref.add_argument('--filter', default=False, help="", required=False, action='store_true')
    parser_pref.add_argument('--surrender', default=False, help="", required=False, action='store_true')
    # -600 for cDNA, pacbio samples
    parser_pref.add_argument('-fp', '--five-prime', type=int, help="", required=False, default=-800)
    parser_pref.add_argument('-tp', '--three-prime', type=int, help="set -1 for deactivation", \
                        required=False, default=-500)
    parser_pref.add_argument('-tcov', '--target-cover', type=int, help="", required=False, default=0.25)

    parser_em = subparsers.add_parser('em', help="arguments for the em module")
    parser_em.add_argument('--pre-init', default=False, \
                        help="use pre-computed abundance / coverage to initalize alpha", \
                        required=False, action='store_true')
    parser_em.add_argument('-e', '--estimate', type=str, help="TSV", \
                        required=False, default=None)
    parser_em.add_argument('-s', '--scores', type=str, help="TSV", \
                        required=False, default=None)
    parser_em.add_argument('-i', '--index', type=str, help="", \
                        required=True)
    parser_em.add_argument('-r', '--rho-thres', type=float, \
                        help="minimum cumulative rho change", \
                        required=False, default=0.0005)
    parser_em.add_argument('-m', '--max-iter', type=int, \
                        help="maximum number of EM iterations", \
                        required=False, default=100)
    parser_em.add_argument('-o', '--out-dir', type=str, help="", \
                        required=True)
    parser_em.add_argument('--drop', default=False, help="", \
                        required=False, action='store_true')
    parser_em.add_argument('--push', default=False, help="", \
                        required=False, action='store_true')
    parser_em.add_argument('-f', '--drop-fac', type=float, \
                        help="factor used to calculate drop threshold", required=False, default=0.3)
    parser_em.add_argument('--naive', default=False, help="", \
                        required=False, action='store_true')
    parser_em.add_argument('-v', '--verbose', default=False, help="", \
                        required=False, action='store_true')
    parser_em.add_argument('--use-score', default=False, help="", \
                        required=False, action='store_true')
    parser_em.add_argument('--relax', default=False, help="", \
                        required=False, action='store_true')
    
    args = parser.parse_args()

    if args.module not in ['align', 'prefilter', 'em']:
        parser.error(f"Invalid module '{args.module}'. Valid options are: align, prefilter, em")
    
    return args

def main():
    args = parse()
    if args.module == "align":
        print("### ALIGN ###")
        align.main(args)
    elif args.module == "prefilter":
        print("### PREFILTER ###")
        prefilter.main(args)
    elif args.module == "em":
        print("### EM ###")
        em.main(args)

if __name__ == "__main__":
    main()