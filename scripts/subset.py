#!/usr/bin/env python

import argparse
import sys
import pickle
import random
import os
from tqdm import tqdm


class Tx:
    def __init__(self, tid, raw):
        self.tid = tid
        self.raw = raw
        # list of exons and CDSs (as gtf lines)
        self.children = list()

    def add_child(self, obj):
        self.children.append(obj)


class Gene:
    def __init__(self, gid, raw):
        self.gid = gid
        self.raw = raw
        # list of Tx objs
        self.children = list()

    def add_child(self, tx):
        self.children.append(tx)


# load gtf file
def load_gtf(fn):
    fh = open(fn, 'r')
    tx_tbl = dict()
    for line in fh:
        if line[0] == '#':
            continue
        clean_ln = line.strip()
        fields = clean_ln.split("\t")
        feature = fields[2]
        infos = fields[8].split(";")
        if feature == "transcript":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tid = val
                    tx = Tx(tid, clean_ln)
                    tx_tbl[tid] = tx
                    break
        elif feature == "exon" or feature == "CDS":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tid = val
                    tx_tbl[tid].add_child(clean_ln)
                    break
    fh.close()
    return tx_tbl


# load gtf w/ genes
def load_gtf_w_genes(fn):
    fh = open(fn, 'r')
    gene_tbl = dict()
    tx_tbl = dict()
    tx2gene_tbl = dict()
    for line in fh:
        if line[0] == '#':
            continue
        clean_ln = line.strip()
        fields = clean_ln.split("\t")
        feature = fields[2]
        infos = fields[8].split(";")
        if feature == "gene":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "gene_id":
                    gid = val
                    gene = Gene(gid, clean_ln)
                    gene_tbl[gid] = gene
                    break
        elif feature == "transcript":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "gene_id":
                    gid = val
                elif key == "transcript_id":
                    tid = val
                    tx = Tx(tid, clean_ln)
                    gene_tbl[gid].add_child(tx)
                    tx_tbl[tid] = tx
                    tx2gene_tbl[tid] = gid
                    break
        elif feature == "exon" or feature == "CDS":
            for info in infos:
                clean_info = info.strip()
                kv_pair = clean_info.split(" ")
                if len(kv_pair) < 2:
                    break
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tid = val
                    tx_tbl[tid].add_child(clean_ln)
                    break
    fh.close()
    return gene_tbl, tx_tbl, tx2gene_tbl


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--index', type=str, help="", required=False, default=None)
    parser.add_argument('-ip', '--index-prefix', type=str, help="", required=False, default=None)
    parser.add_argument('--gtf', type=str, help="", required=True)
    parser.add_argument('--gt-gtf', type=str, help="", required=False, default=None)
    parser.add_argument('--save', default=False, help="", required=False, action='store_true')
    parser.add_argument('-p', type=int, help="", required=True)
    parser.add_argument('-o', '--out-file', type=str, help="", required=True)
    parser.add_argument('--gt-only', default=False, help="", required=False, action='store_true')
    parser.add_argument('--w-genes', default=False, help="", required=False, action='store_true')
    args = parser.parse_args()

    if not 0 <= args.p <= 100:
        print("ARGUMENT ERROR: please provide a percentage argument that is > 1 and <= 100")
        sys.exit()

    if not args.gt_only and args.p == 0:
        print("percentage = 0 is only allowed in gt-only mode")

    if args.w_genes:
        if args.save and args.index_prefix is None:
            print("ARGUMENT ERROR: please provide an index prefix to save the transcriptome")
            sys.exit()
    else:
        if args.save and args.index is None:
            print("ARGUMENT ERROR: please provide an index path to save the transcriptome")
            sys.exit()

    if args.gt_only:
        if args.gt_gtf is None:
            print("ARGUMENT ERROR: please provide the ground truth gtf file to use gt-only mode")
            sys.exit()

    skip = False
    # if index is specified AND it exists, simply load it

    if args.w_genes:
        if args.index_prefix is not None:
            tx_index_fn = os.path.join(args.index_prefix + ".tx.index")
            gene_index_fn = os.path.join(args.index_prefix + ".gene.index")
            tx2gene_index_fn = os.path.join(args.index_prefix + ".tx2gene.index")
            if os.path.exists(tx_index_fn) and os.path.exists(gene_index_fn) and os.path.exists(tx2gene_index_fn):
                with open(tx_index_fn, "rb") as f:
                    tx_tbl = pickle.load(f)
                with open(gene_index_fn, "rb") as f:
                    gene_tbl = pickle.load(f)
                with open(tx2gene_index_fn, "rb") as f:
                    tx2gene_tbl = pickle.load(f)
                skip = True

    else:
        if args.index is not None:
            if os.path.exists(args.index):
                print("loading pre-built transcriptome index...")
                with open(args.index, "rb") as f:
                    tx_tbl = pickle.load(f)
                skip = True
        else:
            print("no pre-built transcriptome index detected... creating it")

    if not skip:
        print("loading the target gtf file...")
        if args.w_genes:
            gene_tbl, tx_tbl, tx2gene_tbl = load_gtf_w_genes(args.gtf)
            if args.save:
                print("saving the target transcriptome as an index...")
                tx_index_fn = os.path.join(args.index_prefix + ".tx.index")
                gene_index_fn = os.path.join(args.index_prefix + ".gene.index")
                tx2gene_index_fn = os.path.join(args.index_prefix + ".tx2gene.index")
                with open(tx_index_fn, "wb") as f:
                    pickle.dump(tx_tbl, f)
                with open(gene_index_fn, "wb") as f:
                    pickle.dump(gene_tbl, f)
                with open(tx2gene_index_fn, "wb") as f:
                    pickle.dump(tx2gene_tbl, f)
        else:
            tx_tbl = load_gtf(args.gtf)
            if args.save:
                print("saving the target transcriptome as an index...")
                with open(args.index, "wb") as f:
                    pickle.dump(tx_tbl, f)

    if args.gt_only:
        gt_tx_tbl = load_gtf(args.gt_gtf)

    if args.gt_only:
        if args.p == 0:
            print("removing all ground truth transcripts")
            del_tid_lst = (gt_tx_tbl.keys())
        else:
            print("picking random indices to remove from the transcriptome")
            gt_tx_num = len(gt_tx_tbl)
            print("total # of ground truth transcripts: " + str(gt_tx_num))
            removal_perc = 100 - args.p
            print("percentage to remove: " + str(removal_perc) + "%")
            removal_num = int(gt_tx_num * removal_perc * 0.01)
            print("# of ground truth transcripts to remove: " + str(removal_num))
            rand_lst = random.sample(range(0, gt_tx_num), removal_num)

            gt_tid_lst = list(gt_tx_tbl.keys())
            del_tid_lst = list()
            for i in tqdm(range(gt_tx_num)):
                tid = gt_tid_lst[i]
                if i in rand_lst:
                    del_tid_lst.append(tid)

    else:
        print("picking random indices to remove from the transcriptome")
        tx_num = len(tx_tbl)
        print("total # of transcripts: " + str(tx_num))
        removal_perc = 100 - args.p
        print("percentage to remove: " + str(removal_perc) + "%")
        removal_num = int(tx_num * removal_perc * 0.01)
        print("# of transcripts to remove: " + str(removal_num))
        rand_lst = random.sample(range(0, tx_num), removal_num)

    print("removing transcripts and writing the kept ones to a gtf file...")
    out_fh = open(args.out_file, 'w')
    if args.w_genes:
        visited_gid_set = set()
    if args.gt_only:
        for tid in tx_tbl:
            if tid not in del_tid_lst:
                if args.w_genes:
                    gid = tx2gene_tbl[tid]
                    gene = gene_tbl[gid]
                    if gid not in visited_gid_set:
                        out_fh.write(gene.raw + "\n")
                        visited_gid_set.add(gid)
                tx = tx_tbl[tid]
                out_fh.write(tx.raw + "\n")
                children = tx.children
                for child in children:
                    out_fh.write(child + "\n")
    else:
        tx_lst = list(tx_tbl.keys())
        for i in tqdm(range(len(tx_lst))):
            if i not in rand_lst:
                tid = tx_lst[i]
                if args.w_genes:
                    gid = tx2gene_tbl[tid]
                    gene = gene_tbl[gid]
                    if gid not in visited_gid_set:
                        out_fh.write(gene.raw + "\n")
                        visited_gid_set.add(gid)
                tx = tx_tbl[tid]
                out_fh.write(tx.raw + "\n")
                children = tx.children
                for child in children:
                    out_fh.write(child + "\n")
    out_fh.close()


if __name__ == "__main__":
    main()
