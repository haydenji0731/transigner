#!/usr/bin/env python

import argparse
import pysam
import os
from time import strftime
import math
import utils


def step_e(assignment):
    for qname in assignment.keys():
        alpha = assignment[qname] # dict() obj
        z = sum(alpha.values())
        for tname in alpha.keys():
            assignment[qname][tname] /= z


def step_m(abundance, assignment):
    rho = dict()
    for tname in abundance.keys():
        rho[tname] = 0
    for qname in assignment:
        alpha = assignment[qname] # dict() object
        for tname in alpha.keys():
            rho[tname] += alpha[tname]
    # normalize
    z = sum(rho.values())
    print("(sanity check) total number of assigned reads: %d" % z)
    for tname in abundance.keys():
        if rho[tname] != 0:
            rho[tname] /= z
    tpm_sum = sum(rho.values())
    print("(sanity check) sum of all rhos: %f" % tpm_sum)
    for qname in assignment:
        alpha = assignment[qname]
        for tname in alpha.keys():
            assignment[tname] = rho[tname]
    return rho


def iter_em(assignment, abundance, min_delta):
    converged = False
    iteration = 0
    while not converged:
        iteration += 1
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "EM iteration #%d:" % iteration)
        step_e(assignment)
        abundance_new = step_m(abundance, assignment)
        converged = has_converged(abundance, abundance_new, min_delta)
        abundance = abundance_new


def has_converged(abundance_old, abundance_new, min_delta):
    max_delta = -math.inf
    for tname in abundance_old.keys():
        delta = abs(abundance_new[tname] - abundance_old[tname])
        max_delta = max(max_delta, delta)
    if max_delta < min_delta:
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Convergence condition satisfied")
        return True
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "delta: %.5f" % max_delta)
    return False


def extract_transcripts(rt):
    abundance = dict()
    fh = open(rt, 'r')
    tot = 0
    for line in fh:
        if line[0] == "#":
            continue
        if line.split("\t")[2] == "transcript":
            tot += 1
            fields = line.split("\t")[8].split(";")
            for field in fields:
                k = field.strip().split(" ")[0]
                v = field.strip().split(" ")[1].replace('"', '')
                if k == "transcript_id":
                    abundance[v] = 0
                    break
    for tname in abundance.keys():
        abundance[tname] = 1 / tot
    return abundance, tot


def assign_compatibility(in_aln, tot_trs):
    assignment = dict()
    unmapped = 0
    qnames = set()
    with pysam.AlignmentFile(in_aln, 'rb') as fh:
        for brec in fh:
            if brec.is_unmapped:
                unmapped += 1
                continue
            else:
                qname = brec.query_name
                qnames.add(qname)
    mapped = len(qnames)
    for name in qnames:
        assignment[name] = dict()
    with pysam.AlignmentFile(in_aln, 'rb') as fh:
        for brec in fh:
            if brec.is_unmapped:
                continue
            else:
                qname = brec.query_name
                tname = brec.reference_name
                assignment[qname][tname] = 1 / tot_trs
    tot = unmapped + mapped
    print("Loaded total %d reads of which %d are mapped / %d are unmapped" % (tot, mapped, unmapped))
    return assignment, mapped


def main():
    parser = argparse.ArgumentParser(description="run EM algorithm to (1) quantify transcript abundance and (2) "
                                                 "estimate read membership likelihoods")
    parser.add_argument('-i', '--input_aln', type=str, help="input alignment file", required=True)
    parser.add_argument('-ref-gtf', '--ref_gtf', type=str, help="reference transcriptome annotation to match against",
                        required=True)
    parser.add_argument('-thres', '--threshold', type=int, help="threshold for stopping em", default=0.001,
                        nargs='?', const=0.001)
    parser.add_argument('-o', '--output_dir', type=str, help="output directory", required=True)
    parser.add_argument('-op', '--output_prefix', type=str, help="output files prefix", default='quant',
                        nargs='?', const='quant')

    args = parser.parse_args()
    rt = args.ref_gtf
    in_aln = args.input_aln
    min_delta = args.threshold
    out_dir = args.output_dir
    out_prefix = args.output_prefix

    abundance, tot_trs = extract_transcripts(rt)

    if not os.path.exists(in_aln + ".bai"):
        pysam.index(in_aln)

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading query alignment file")
    assignment, tot_mapped = assign_compatibility(in_aln, tot_trs)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished loading. Beginning EM.")
    iter_em(assignment, abundance, min_delta)
    utils.write_results(out_dir, out_prefix, assignment, abundance, tot_mapped)


if __name__ == "__main__":
    main()
