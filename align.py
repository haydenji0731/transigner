#!/usr/bin/env python

import argparse
from subprocess import call
from time import strftime


def align(in_fasta, rt, threads, op, odir):
    # align as if genomic reads
    minimap_cmd = "minimap2 -ax map-ont -t " + str(threads) + " " + rt + " " + in_fasta + \
                  " | samtools sort -o " + odir + "/" + op + ".bam"
    print(minimap_cmd)
    call(minimap_cmd, shell=True)
    samtools_cmd = "samtools index " + odir + "/" + op + ".bam -@ " + str(threads)
    print(samtools_cmd)
    call(samtools_cmd, shell=True)


def main():
    
    parser = argparse.ArgumentParser(description="align reads to the transcriptome")
    parser.add_argument('-i', '--input_fasta', type=str, help="input fasta file containing reads", required=True)
    parser.add_argument('-ref-fa', '--ref_fasta', type=str, help="reference transcriptome fasta to align to",
                        required=True)
    parser.add_argument('-o', '--output_dir', type=str, help="alignment output directory", required=True)
    parser.add_argument('-op', '--output_prefix', type=str, help="output alignment file prefix", default='aln',
                        nargs='?', const='aln')
    parser.add_argument('-t', '--threads', type=int, help="number of threads used in alignment", default=1, nargs='?',
                        const=1)
    args = parser.parse_args()
    in_fasta = args.input_fasta
    rt = args.ref_fasta
    threads = args.threads
    odir = args.output_dir
    op = args.output_prefix
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Starting alignment to the reference transcriptome using minimap2")
    align(in_fasta, rt, threads, op, odir)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Alignment completed.")


if __name__ == "__main__":
    main()
