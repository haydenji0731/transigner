#!/usr/bin/env python

import argparse
from subprocess import call
from time import strftime


def align(in_fasta, ref_t, aln_fn, threads, of):
    if of == "bam":
        minimap_cmd = "minimap2 -ax map-ont -t " + str(threads) + " " + ref_t + " " + in_fasta + \
                      " | samtools view -b > " + aln_fn
        print(minimap_cmd)
        # call(minimap_cmd)
    else:
        minimap_cmd = "minimap2 -ax map-ont -t " + str(threads) + " " + ref_t + " " + in_fasta + " > " + aln_fn
        print(minimap_cmd)
        # call(minimap_cmd)


def main():
    parser = argparse.ArgumentParser(description="align reads to the transcriptome")
    parser.add_argument('-i', '--input_fasta', type=str, help="input fasta file containing reads", required=True)
    parser.add_argument('-rt', '--ref_transcriptome', type=str, help="reference transcriptome to align reads to",
                        required=True)
    parser.add_argument('-o', '--output_dir', type=str, help="alignment output directory", required=True)
    parser.add_argument('-op', '--output_prefix', type=str, help="output alignment file prefix", default='aln',
                        nargs='?', const='aln')
    parser.add_argument('-of', '--output_format', type=str, help="output alignment file format (bam/sam)",
                        default='bam', nargs='?', const='cram')
    parser.add_argument('-t', '--threads', type=int, help="number of threads used in alignment", default=1, nargs='?',
                        const=1)
    args = parser.parse_args()
    if args.output_format.lower() == "bam":
        of = 'bam'
        aln_fn = args.output_prefix + ".bam"
    else:
        of = 'sam'
        aln_fn = args.output_prefix + ".sam"

    in_fasta = args.input_fasta
    ref_t = args.ref_transcriptome
    threads = args.threads
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Starting alignment to the reference transcriptome using minimap2")
    align(in_fasta, ref_t, aln_fn, threads, of)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Alignment completed.")


if __name__ == "__main__":
    main()
