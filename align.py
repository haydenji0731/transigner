#!/usr/bin/env python

import argparse

# def align(in_fasta, num_threads, prefix):
#     aln_fn = prefix + "_aln.bam"
#     call("")


def main():
    parser = argparse.ArgumentParser(description="align reads to the transcriptome")
    parser.add_argument('-i', '--input_fasta', type=str, help="input fasta file containing reads", required=True)
    parser.add_argument('-rt', '--ref_transcriptome', type=str, help="reference transcriptome to align reads to", required=True)
    parser.add_argument('-o', '--output_dir', type=str, help="alignment output directory", required=True)
    parser.add_argument('-of', '--output_prefix', type=str, help="output alignment file prefix", default='aln')
    args = parser.parse_args()


if __name__ == "__main__":
    main()