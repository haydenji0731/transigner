#!/usr/bin/env python

import csv
import math


def write_results(out_dir, out_prefix, assignment, abundance, tot_mapped):
    abundance_fn = out_dir + "/" + out_prefix + ".abundance.tsv"
    assignment_fn = out_dir + "/" + out_prefix + ".assignment.tsv"
    with open(abundance_fn, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['transcript_id', 'relative abundance', 'raw read count', 'TPM'])
        for tname in abundance.keys():
            rho = abundance[tname]
            rc = rho * tot_mapped
            tpm = rho * 1000000
            writer.writerow([tname, rho, rc, tpm])

    with open(assignment_fn, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        # TODO: check confidence
        writer.writerow(['read_id', 'transcript_id', 'confidence'])
        for rname in assignment.keys():
            alpha = assignment[rname]
            max_l = -math.inf
            for tname in alpha.keys():
                if alpha[tname] > max_l:
                    max_l = alpha[tname]
                    best_match = (tname, max_l)
            writer.writerow([rname, best_match[0], best_match[1]])


