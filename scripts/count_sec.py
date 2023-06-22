#!/usr/bin/env python

import pysam
import sys
from time import strftime


def count(aln_fn):
    # initialize
    aln_d = dict()
    unmapped = 0
    qname_prev = ""

    with pysam.AlignmentFile(aln_fn, 'rb') as fh:
        for brec in fh:
            if not brec.is_unmapped:
                qname = brec.query_name
                if qname != qname_prev:
                    # previous read is NIL
                    if qname_prev != "":
                        aln_d[qname_prev] = tot_sec
                    tot_sec = 0
                    qname_prev = qname
                if brec.is_secondary:
                    tot_sec += 1
            else:
                unmapped += 1
    print("total # of unmapped reads: " + str(unmapped))
    return aln_d


def compare(d1, d2):
    identical = True
    for query in d1.keys():
        d1_val = d1[query]
        try:
            d2_val = d2[query]
        except:
            print("WARNING: a mismatch in query set detected")
        if d1_val != d2_val:
            identical = False
            print(query)
            if d1_val < d2_val:
                print("alignment #1 is smaller")
            else:
                print("alignment #2 is smaller")
            break
    return identical


def main(aln_fn_1, aln_fn_2):
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Processing alignment file: " + aln_fn_1)
    aln_d_1 = count(aln_fn_1)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Processing alignment file: " + aln_fn_2)
    aln_d_2 = count(aln_fn_2)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Comparing 2 alignment files")
    is_identical = compare(aln_d_1, aln_d_2)
    print("These 2 alignment files have same # of secondary alignments: " + str(is_identical))
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Run completed")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])