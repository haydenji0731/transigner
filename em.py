#!/usr/bin/env python

import argparse
import pysam
import os
from time import strftime
import utils
import math


assignment = dict()
abundance = dict()


def step_e():
    global assignment
    for qname in assignment.keys():
        alpha = assignment[qname] # dict() obj
        z = sum(alpha.values())
        for tname in alpha.keys():
            assignment[qname][tname] /= z


def step_m():
    global abundance
    global assignment
    rho = dict()
    for tname in abundance.keys():
        rho[tname] = 0
    for qname in assignment.keys():
        alpha = assignment[qname] # dict() object
        for tname in alpha.keys():
            rho[tname] += alpha[tname]
    z = sum(rho.values())
    # print("(sanity check) total number of assigned reads: %d" % z)
    for tname in abundance.keys():
        if rho[tname] != 0:
            rho[tname] /= z
    # rho_sum = sum(rho.values())
    # print("(sanity check) sum of all rhos: %f" % rho_sum)
    return rho


def iter_em(min_delta, max_iter):
    converged = False
    iteration = 0
    global abundance
    global assignment
    while not converged:
        iteration += 1
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "EM iteration #%d:" % iteration)
        for qname in assignment.keys():
            alpha = assignment[qname]
            for tname in alpha.keys():
                assignment[qname][tname] = abundance[tname]
        step_e()
        abundance_new = step_m()
        converged = has_converged(abundance, abundance_new, min_delta)
        abundance = abundance_new
        if iteration >= max_iter:
            converged = True


def has_converged(abundance_old, abundance_new, min_delta):
    max_delta = -math.inf
    for tname in abundance_old.keys():
        delta = abs(abundance_new[tname] - abundance_old[tname])
        max_delta = max(max_delta, delta)
    max_delta *= 1000000
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "delta: %.10f" % max_delta)
    if max_delta < min_delta:
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Convergence condition satisfied")
        return True
    return False


def extract_transcripts(rt):
    global abundance
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


def assign_compatibility(in_aln):
    global assignment
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
                assignment[qname][tname] = 0
    tot = unmapped + mapped
    print("Loaded total %d reads of which %d are mapped / %d are unmapped" % (tot, mapped, unmapped))
    return mapped


def main():
    parser = argparse.ArgumentParser(description="run EM algorithm to (1) quantify transcript abundance and (2) "
                                                 "estimate read membership likelihoods")
    parser.add_argument('-i', '--input_aln', type=str, help="input alignment file", required=True)
    parser.add_argument('-ref-gtf', '--ref_gtf', type=str, help="reference transcriptome annotation to match against",
                        required=True)
    parser.add_argument('-thres', '--threshold', type=float, help="min TPM change for stopping EM", default=10)
    parser.add_argument('-max-iter', '--max_iteration', type=int, help="maximum number of EM iterations", default=100)
    parser.add_argument('-o', '--output_dir', type=str, help="output directory", required=True)
    parser.add_argument('-op', '--output_prefix', type=str, help="output files prefix", default='quant')

    args = parser.parse_args()
    rt = args.ref_gtf
    in_aln = args.input_aln
    min_delta = args.threshold
    out_dir = args.output_dir
    out_prefix = args.output_prefix
    max_iter = args.max_iteration

    print("\n-----------------------")
    print("Run Parameters")
    print("-----------------------")
    print("reference gtf: " + rt)
    print("input bam: " + in_aln)
    print("minimum TPM change: " + str(min_delta))
    print("maximum iteration: " + str(max_iter))
    print("output directory: " + out_dir)
    print("output prefix: " + out_prefix)
    print("-----------------------")

    extract_transcripts(rt)

    if not os.path.exists(in_aln + ".bai"):
        pysam.index(in_aln)

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading query alignment file")
    tot_mapped = assign_compatibility(in_aln)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished loading. Beginning EM.")
    iter_em(min_delta, max_iter)
    utils.write_results(out_dir, out_prefix, assignment, abundance, tot_mapped)


if __name__ == "__main__":
    main()
