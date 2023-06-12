#!/usr/bin/env python

import argparse
from subprocess import call
from time import strftime

g2t_d = dict()


def align(in_fastq, rt, threads, op, odir, sec_num):
    # align as if genomic reads
    minimap_cmd = "minimap2 -ax map-ont -N " + str(sec_num) + " -t " + str(threads) + " " + rt + " " + in_fastq + \
                  " | samtools sort -o " + odir + "/" + op + ".bam"
    print(minimap_cmd)
    call(minimap_cmd, shell=True)
    samtools_cmd = "samtools index " + odir + "/" + op + ".bam -@ " + str(threads)
    print(samtools_cmd)
    call(samtools_cmd, shell=True)


def get_tx_per_gene(ref_gtf):
    global g2t_d
    max_tpg = 0
    ref_gtfh = open(ref_gtf, 'r')
    for line in ref_gtfh:
        if line[0] == '#':
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            infos = fields[8].split(';')
            for info in infos:
                if info == '\n':
                    break
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "ref_gene_id":
                    gid = val
                    if gid not in g2t_d.keys():
                        g2t_d[gid] = 0
                    g2t_d[gid] += 1
                    break

        # if feature == "gene":
        #     infos = fields[8].split(';')
        #     for info in infos:
        #         info_clean = info.strip()
        #         kv_pair = info_clean.split(" ")
        #         key = kv_pair[0]
        #         val = kv_pair[1].replace('"', '')
        #         if key == "ref_gene_id":
        #             gid = val
        #             g2t_d[gid] = 0
        #             break
        # elif feature == "transcript":
        #     g2t_d[gid] += 1
    for gene in g2t_d.keys():
        num_tx = g2t_d[gene]
        if max_tpg < num_tx:
            max_gid = gene
        max_tpg = max(max_tpg, num_tx)
    return max_gid, max_tpg


def main():
    
    parser = argparse.ArgumentParser(description="align reads to the transcriptome")
    parser.add_argument('-i', '--input_fastq', type=str, help="input fastq file containing reads", required=True)
    parser.add_argument('-ref-fa', '--ref_fasta', type=str, help="reference transcriptome fasta to align to",
                        required=True)
    parser.add_argument('--ref-gtf', type=str, help="reference transcriptome gtf", required=False, default="")
    parser.add_argument('-o', '--output_dir', type=str, help="alignment output directory", required=True)
    parser.add_argument('-op', '--output_prefix', type=str, help="output alignment file prefix", default='aln',
                        required=False)
    parser.add_argument('-sN', '--secondary_num', type=int, help="maximum number of secondary alignments", default=10,
                        required=False)
    parser.add_argument('-t', '--threads', type=int, help="number of threads used in alignment", default=1,
                        required=False)
    args = parser.parse_args()
    in_fastq = args.input_fastq
    rt = args.ref_fasta
    threads = args.threads
    odir = args.output_dir
    op = args.output_prefix
    if args.ref_gtf != "":
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Reference annotation was provided. Extracting the maximum number of "
                                                 "transcripts per gene locus.")
        gid, tpg = get_tx_per_gene(args.ref_gtf)
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Maximum number of transcripts per gene locus: " + str(tpg))
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "The corresponding gene locus: " + gid)
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Maximum number of secondary alignments will be ignored.")
        # TODO: try adding some padding
        sec_num = tpg + 100
        #sec_num = tpg
    else:
        sec_num = args.secondary_num
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Starting alignment to the reference transcriptome using minimap2")
    align(in_fastq, rt, threads, op, odir, sec_num)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Alignment completed.")


if __name__ == "__main__":
    main()
