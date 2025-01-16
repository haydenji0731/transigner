#!/usr/bin/env python

from transigner.utils import *
from transigner import align, pre, em

def parse():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--version', action='version', version='%(prog)s v0.0.1')
    subparsers = parser.add_subparsers(dest='module', \
                               help="[align, pre, em]")
    # run all modules; easy peasy
    # align module        
    parser_align = subparsers.add_parser('align', help="")
    parser_align.add_argument('-q', '--query', type=str, help="", required=True)
    parser_align.add_argument('-t', '--target', type=str, help="", required=True)
    parser_align.add_argument('-d', '--out-dir', type=str, help="", required=False, default=".")
    parser_align.add_argument('-o', '--out-file', type=str, help="", required=False, default="temp.bam")
    parser_align.add_argument('-n', type=int, help="", default=181, required=False)
    parser_align.add_argument('-p', '--threads', type=int, help="", default=1, \
                        required=False)
    parser_align.add_argument('--preset', type=str, help="", \
                            required=False, default="map-ont")
    parser_align.add_argument('--dev', default=False, help="", \
                            required=False, action='store_true')
    parser_align.add_argument('--index-opts', type=str, help="", \
                            required=False, default="")
    parser_align.add_argument('--map-opts', type=str, help="", \
                            required=False, default="")
    parser_align.add_argument('--base-aln-opts', type=str, help="", \
                            required=False, default="")
    # pre module
    parser_pre = subparsers.add_parser('pre', help="")
    parser_pre.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser_pre.add_argument('-d', '--out-dir', type=str, help="", required=True)
    parser_pre.add_argument('-e', '--estimates', type=str, help="", required=False, default=None)
    parser_pre.add_argument('--dp-score-model', nargs=1, type=str, help="", \
                            required=False, default='e', choices=['e', 'pl'])
    parser_pre.add_argument('--dsm-opts', nargs=2, type=float, help="amplitude, decay rate",\
                            required=False, default=[0.0, 5.0]) # default 10.0
    parser_pre.add_argument('--use-filter', default=False, help="", required=False, action='store_true')
    parser_pre.add_argument('--filt-opts', nargs=3, type=float, help="5' 3' tcov", \
                            required=False, default=[0, 0, 0])
    parser_pre.add_argument('--use-psw', default=False, help="", required=False, action='store_true')
    parser_pre.add_argument('--dev', default=False, help="", required=False, action='store_true')
    parser_pre.add_argument('--spiked', default=False, help="", required=False, action='store_true')

    # em module
    parser_em = subparsers.add_parser('em', help="")
    parser_em.add_argument('-s', '--scores', type=str, help="", required=True)
    parser_em.add_argument('-d', '--out-dir', type=str, help="", required=True)
    parser_em.add_argument('-u', '--unmapped', type=str, help="", required=True)
    parser_em.add_argument('-m', '--tmap-file', type=str, help="", required=True)
    parser_em.add_argument('-n', '--num-iter', type=int, help="", required=False, default=1000)
    parser_em.add_argument('-c', '--cvrg-thres', help="minimum delt", \
                           required=False, default='auto')
    parser_em.add_argument('-utol', '--unmapped-tol', nargs=2, type=float, help="", \
                           required=False, default=[3.5, 3.9])
    parser_em.add_argument('-dtype', type=str, help="[ont_drna, ont_cdna, pacbio]", \
                           required=False, default='ont_drna')
    parser_em.add_argument('-r', '--relax-thres', type=float, help="", required=False, default=0.1)
    parser_em.add_argument('-p', '--threads', type=int, help="", required=False, default=1)
    parser_em.add_argument('--naive', default=False, help="", required=False, \
                           action='store_true')
    parser_em.add_argument('--dev', default=False, help="", required=False, \
                           action='store_true')
    parser_em.add_argument('--push', default=False, help="", required=False, \
                           action='store_true')
    # TODO: make this default
    parser_em.add_argument('--drop', default=False, help="", required=False, \
                           action='store_true')
    parser_em.add_argument('--relax', default=False, help="", required=False, \
                           action='store_true')
    parser_em.add_argument('-df', '--drop-fac', type=float, help="", \
                           required=False, default=0.1)

    args = parser.parse_args()
    if args.module not in ['align', 'pre', 'em']:
        parser.error(f"Invalid module '{args.module}'. Valid options are: align, pre, em")
    return args

def main() -> None:
    args = parse()
    if args.module == "align":
        align.main(args)
    elif args.module == "pre":
        pre.main(args)
    elif args.module == "em":
        em.main(args)
