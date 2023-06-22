#!/usr/bin/env python

import sys
import pyfastx


def get_gt_set(reads_fn):
    tx_set = set()
    fq = pyfastx.Fastx(reads_fn)
    for name, seq, _ in fq:
        splits = name.split("_")
        if splits[0] == "unassigned":
            gt = splits[0] + "_" + splits[1] + "_" + splits[2]
        else:
            gt = splits[0] + "_" + splits[1]
        tx_set.add(gt)
    return tx_set


def match_tx2gene(gtf_fn, tx_set):
    gene_set = set()
    gtf_fh = open(gtf_fn)
    for line in gtf_fh:
        if line[0] == '#':
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            infos = fields[8].split(";")
            for info in infos:
                info = info.strip()
                kv_pair = info.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "gene_id":
                    gid = val
                elif key == "transcript_id":
                    tid = val.split(".")[0]
                    break
            if tid in tx_set:
                gene_set.add(gid)
    gtf_fh.close()
    return gene_set


def filter_gtf(gtf_fn, out_prefix, tx_set, gene_set):
    gtf_fh = open(gtf_fn, 'r')
    kept_fn = out_prefix + "_kept.gtf"
    kept_fh = open(kept_fn, 'w')
    removed_fn = out_prefix + "_removed.gtf"
    removed_fh = open(removed_fn, 'w')
    if gene_set is None:
        for line in gtf_fh:
            if line[0] == '#':
                continue
            fields = line.split("\t")
            infos = fields[8].split(";")
            for info in infos:
                info = info.strip()
                kv_pair = info.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tx_id = val.split(".")[0]
                    if tx_id in tx_set:
                        removed_fh.write(line)
                    else:
                        kept_fh.write(line)
                    break
    else:
        for line in gtf_fh:
            if line[0] == '#':
                continue
            fields = line.split("\t")
            infos = fields[8].split(";")
            for info in infos:
                info = info.strip()
                kv_pair = info.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "gene_id":
                    gene_id = val
                    if gene_id in gene_set:
                        removed_fh.write(line)
                    else:
                        kept_fh.write(line)
    gtf_fh.close()
    removed_fh.close()
    kept_fh.close()


def main(gtf_fn, reads_fn, out_prefix, with_genes):
    tx_set = get_gt_set(reads_fn)
    if with_genes:
        gene_set = match_tx2gene(gtf_fn, tx_set)
    else:
        gene_set = None
    filter_gtf(gtf_fn, out_prefix, tx_set, gene_set)


if __name__ == "__main__":
    # gtf_fn, reads_fn, out_prefix, chrom, with_genes
    main(sys.argv[1], sys.argv[2], sys.argv[3], bool(sys.argv[4]))
