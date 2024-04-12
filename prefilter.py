#!/usr/bin/env python

import pysam
from tqdm import tqdm
from subprocess import call
import argparse
import pyfastx
import os
import pickle
import sys
from utils import build_tbls_paf, build_score_tbl
import json


def load_cmpt_tbl(cmpt_fn, reads_num, tx_tbl):
    cmpt_fh = open(cmpt_fn, 'r')
    cmpt_tbl = [set() for _ in range(reads_num)]
    ctr = -1
    reads_lst = list()
    for line in cmpt_fh:
        fields = line.strip().split("\t")
        for i in range(len(fields)):
            field = fields[i]
            if i == 0:
                ctr += 1
                qname = field
                reads_lst.append(qname)
            else:
                tname = field
                t_idx = tx_tbl[tname]
                cmpt_tbl[ctr].add(t_idx)
    cmpt_fh.close()
    return cmpt_tbl, reads_lst


def load_score_tbl(score_fn, reads_num, tx_tbl):
    score_tbl = [dict() for _ in range(reads_num)]
    score_fh = open(score_fn, 'r')
    ctr = -1
    reads_lst = list()
    for line in score_fh:
        fields = line.strip().split("\t")
        for i in range(len(fields)):
            field = fields[i]
            if i == 0:
                ctr += 1
                qname = field
                reads_lst.append(qname)
            else:
                tname = field.split(",")[0].replace('(', '')
                score = float(field.split(",")[1].strip().replace(')', ''))
                ti = tx_tbl[tname]
                score_tbl[ctr][ti] = score
    score_fh.close()
    return score_tbl, reads_lst


def load_transcriptome(gtf_fn):
    tx_lst = list()
    gtf_fh = open(gtf_fn, 'r')
    for line in gtf_fh:
        if line[0] == '#':
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            infos = fields[8].split(";")
            for info in infos:
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                if kv_pair[0] == "transcript_id":
                    tid = kv_pair[1].replace('"', '')
                    break
            tx_lst.append(tid)
    gtf_fh.close()
    return tx_lst


def subset_target(cmpt_tbl, tx_lst, reads_lst):
    # TODO:
    reads_num = len(reads_lst)
    cmpt_tx_set = set()
    for i in tqdm(range(reads_num)):
        t_idx_set = cmpt_tbl[i]
        for t_idx in t_idx_set:
            tname = tx_lst[t_idx]
            cmpt_tx_set.add(tname)
    return cmpt_tx_set


def compute_tlen(target_fa):
    tlen_tbl = dict()
    for seq in tqdm(target_fa):
        if seq.name in tlen_tbl:
            print("how can this happen?")
        tlen_tbl[seq.name] = len(seq)
    return tlen_tbl


def load_aln(aln_fn, tlen_tbl):
    reads_tbl = dict()
    primary_tbl = dict()
    primary_st_tbl = dict()
    primary_en_tbl = dict()
    with pysam.AlignmentFile(aln_fn, 'rb') as fh:
        for brec in tqdm(fh):
            qname = brec.query_name
            if brec.is_unmapped:
                continue
            else:
                if qname in reads_tbl:
                    reads_tbl[qname].append(brec)
                else:
                    reads_tbl[qname] = [brec]
                if not brec.is_secondary and not brec.is_supplementary:
                    score = int(brec.get_tag("ms"))
                    primary_tbl[qname] = (brec.reference_name, score)
                    primary_st_tbl[qname] = (brec.reference_name, brec.reference_start)
                    # TODO: does this work? :0
                    tlen = tlen_tbl[brec.reference_name]
                    # record the distance from the 3' end
                    primary_en_tbl[qname] = (brec.reference_name, tlen - brec.reference_end)
    return reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl


def build_cmpt_tbl_bam(reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl, three_thres, five_thres,
                       tx_tbl, tlen_tbl, reads_lst, filter_by_pos):
    reads_num = len(reads_tbl)

    cmpt_tbl = [set() for _ in range(reads_num)]
    ratio_tbl = [dict() for _ in range(reads_num)]

    for i in tqdm(range(reads_num)):
        qname = reads_lst[i]
        brecs = reads_tbl[qname]
        for brec in brecs:

            tname = brec.reference_name

            if filter_by_pos:
                pri_st_pos = primary_st_tbl[qname][1]
                # skip if the 5' mapping position too far to the right of the primary
                if (pri_st_pos - brec.reference_start) < five_thres:
                    continue

                # TODO: remove this stub as well as other dependencies
                #pri_en_dist = primary_en_tbl[qname][1]
                #tlen = tlen_tbl[brec.reference_name]

                #if (pri_en_dist - (tlen - brec.reference_end)) < three_thres:
                #    continue

            t_idx = tx_tbl[tname]
            # add to the cmpt set for this read
            cmpt_tbl[i].add(t_idx)

            score = int(brec.get_tag("ms"))
            pri_score = primary_tbl[qname][1]
            ratio = score / pri_score

            if t_idx in ratio_tbl[i]:
                ratio_tbl[i][t_idx] = max(ratio_tbl[i][t_idx], ratio)
            else:
                ratio_tbl[i][t_idx] = ratio

    return cmpt_tbl, ratio_tbl


def build_cmpt_tbl_bam_decay(reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl, gamma_3, beta_3, gamma_5, beta_5,
                             tx_tbl, tlen_tbl, reads_lst):
    reads_num = len(reads_tbl)

    cmpt_tbl = [set() for _ in range(reads_num)]
    ratio_tbl = [dict() for _ in range(reads_num)]

    for i in tqdm(range(reads_num)):
        qname = reads_lst[i]
        brecs = reads_tbl[qname]
        for brec in brecs:

            tname = brec.reference_name
            pri_st_pos = primary_st_tbl[qname][1]
            rel_d = max(pri_st_pos - brec.reference_start, 0)

            five_bias = (1 - gamma_5) ** (max(rel_d - beta_5, 0))

            # pri_en_pos = primary_en_tbl[qname][1]
            tlen = tlen_tbl[brec.reference_name]
            d_from_three = tlen - brec.reference_end
            pri_en_dist = primary_en_tbl[qname][1]
            rel_d = max(pri_en_dist - d_from_three, 0)

            three_bias = (1 - gamma_3) ** (max(rel_d - beta_3, 0))

            t_idx = tx_tbl[tname]
            # add to the cmpt set for this read
            cmpt_tbl[i].add(t_idx)

            score = int(brec.get_tag("ms"))
            pri_score = primary_tbl[qname][1]
            ratio = score / pri_score * five_bias * three_bias

            if t_idx in ratio_tbl[i]:
                ratio_tbl[i][t_idx] = max(ratio_tbl[i][t_idx], ratio)
            else:
                ratio_tbl[i][t_idx] = ratio

    return cmpt_tbl, ratio_tbl


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-query', '--query-fastq', type=str, help="", required=True)
    parser.add_argument('-target', '--target-fasta', type=str, help="", required=True)
    parser.add_argument('-aln', type=str, help="", required=True)
    parser.add_argument('-gtf', type=str, help="", required=True)
    # parser.add_argument('-jf', '--jf-aligner', type=str, help="", required=False, default=None)
    # TODO: test this parameter
    parser.add_argument('--align', default=False, help="", required=False, action='store_true')
    parser.add_argument('-tmp', '--tmp-dir', type=str, help="", required=True)
    parser.add_argument('-t', '--threads', type=int, help="", required=False, default=1)
    parser.add_argument('-k', '--kmer-size', type=int, help="", required=False, default=20)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=True)
    parser.add_argument('--load', default=False, help="", required=False, action='store_true')
    parser.add_argument('-i', '--index-path', type=str, help="", required=False, default=None)
    parser.add_argument('--score-ratio', type=float, help="", required=False, nargs="+")
    parser.add_argument('--five-prime', type=float, help="", required=False, default=-800)
    parser.add_argument('--three-prime', type=float, help="", required=False, default=-150)
    parser.add_argument('--filter', default=False, help="", required=False, action='store_true')
    parser.add_argument('--decay', default=False, help="", required=False, action='store_true')

    parser.add_argument('-g3', '--gamma-3', type=float, help="", required=False, default=0.01)
    parser.add_argument('-b3', '--beta-3', type=float, help="", required=False, default=100)
    parser.add_argument('-g5', '--gamma-5', type=float, help="", required=False, default=0.01)
    parser.add_argument('-b5', '--beta-5', type=float, help="", required=False, default=600)
    parser.add_argument('--keep', default=False, help="", required=False, action='store_true')

    parser.add_argument('--skip-jf', default=False, help="", required=False, action='store_true')
    args = parser.parse_args()

    args_fn = os.path.join(args.out_dir, "pref_cmd.txt")
    with open(args_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    print("loading target transcriptome")
    target_fa = pyfastx.Fasta(args.target_fasta)

    print("extracting target transcript lengths")
    tlen_tbl = compute_tlen(target_fa)

    # print jf_aligner warning
    if not args.skip_jf:
        print("jf_aligner option activated.")
        print("please make sure you have jf_aligner path in your PATH variable")

    # TODO: fix this
    if args.load:
        if args.index_path is None:
            print("no transcriptome index path provided. exiting.")
            sys.exit()
        tx_index = os.path.join(args.index_path, "transcriptome.index")
        print("Loading target transcriptome at", tx_index)
        with open(tx_index, "rb") as f:
            tx_tbl = pickle.load(f)
        tmp_lst = list()
        for tname in tx_tbl:
            tmp_lst.append((tname, tx_tbl[tname]))
        sorted_tmp_lst = sorted(tmp_lst, key=lambda x: x[1])
        tx_lst = [x for (x, i) in sorted_tmp_lst]
        print("Loading mm2 compatibility matrix and score tables")
        mm2_cmpt_fn = os.path.join(args.index_path, "cmpt_tbl.tsv")
        cmpt_tbl, reads_lst = load_cmpt_tbl(mm2_cmpt_fn, 8920552, tx_tbl)
        mm2_ms_fn = os.path.join(args.index_path, "ms_tbl.tsv")
        ms_tbl, reads_lst = load_score_tbl(mm2_ms_fn, 8920552, tx_tbl)
        print("Loading jf compatibility matrix and score tables")
        jf_cmpt_fn = os.path.join(args.index_path, "cmpt_tbl_jf.tsv")
        jf_cover_fn = os.path.join(args.index_path, "cover_tbl.tsv")
        cmpt_tbl_jf, reads_lst_jf = load_cmpt_tbl(jf_cmpt_fn, 8252778, tx_tbl)
        cover_tbl, reads_lst_jf = load_score_tbl(jf_cover_fn, 8252778, tx_tbl)

    else:
        print("loading alignment")
        reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl = load_aln(args.aln, tlen_tbl)
        reads_lst = list(reads_tbl.keys())
        print("loading target transcriptome annotation")
        tx_lst = load_transcriptome(args.gtf)

        tx_tbl = dict()
        for i in range(len(tx_lst)):
            tname = tx_lst[i]
            tx_tbl[tname] = i

        tx_index = os.path.join(args.out_dir, "transcriptome.index")
        print("Saving transcriptome index to", tx_index)
        with open(tx_index, "wb") as f:
            pickle.dump(tx_tbl, f)

        print("building compatibility matrix")
        if args.decay:
            cmpt_tbl, ms_tbl = build_cmpt_tbl_bam_decay(reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl,
                                                        args.gamma_3, args.beta_3, args.gamma_5, args.beta_5, tx_tbl,
                                                        tlen_tbl, reads_lst)
        else:
            cmpt_tbl, ms_tbl = build_cmpt_tbl_bam(reads_tbl, primary_tbl, primary_st_tbl, primary_en_tbl,
                                                  args.three_prime,
                                                  args.five_prime, tx_tbl, tlen_tbl, reads_lst, args.filter)

        print("writing out compatibility matrix")
        # print out the compatibility matrix (excluding ratios)
        cmpt_fn = os.path.join(args.out_dir, "cmpt_tbl.tsv")
        cmpt_fh = open(cmpt_fn, 'w')
        reads_num = len(reads_lst)
        for i in tqdm(range(reads_num)):
            qname = reads_lst[i]
            cmpt_fh.write(qname)
            t_idx_set = cmpt_tbl[i]
            for t_idx in t_idx_set:
                tname = tx_lst[t_idx]
                cmpt_fh.write("\t" + tname)
            cmpt_fh.write("\n")
        cmpt_fh.close()

        if os.path.exists(args.tmp_dir):
            print("tmp directory exists. removing it to start clean.")
            rm_cmd = "rm -rf " + args.tmp_dir
            call(rm_cmd, shell=True)

        if args.skip_jf:
            print("writing out the score matrix")
            scores_fn = os.path.join(args.out_dir, "scores.tsv")
            scores_fh = open(scores_fn, 'w')
            for i in tqdm(range(len(reads_lst))):
                qname = reads_lst[i]
                scores_fh.write(qname)
                ti_set = cmpt_tbl[i]
                for ti in ti_set:
                    score = ms_tbl[i][ti]
                    tname = tx_lst[ti]
                    scores_fh.write("\t(" + tname + ", " + str(score) + ")")
                scores_fh.write("\n")
            scores_fh.close()
        else:
            mkdir_cmd = "mkdir " + args.tmp_dir
            call(mkdir_cmd, shell=True)

            print("subsetting target transcriptome")
            cmpt_tx_set = subset_target(cmpt_tbl, tx_lst, reads_lst)
            subset_target_fn = os.path.join(args.tmp_dir, "cmpt_txs.fa")
            subset_target_fh = open(subset_target_fn, 'w')
            for tname in cmpt_tx_set:
                tx = target_fa[tname]
                subset_target_fh.write(tx.raw)
            subset_target_fh.close()

            aln_fn = os.path.join(args.out_dir, "jf_aln.txt")
            aln_cmd = "cat " + args.query_fastq + " | jf_aligner -t " + str(args.threads) + " -B 17.0 -H -s " + \
                      str(len(cmpt_tx_set)) + " -m " + str(
                args.kmer_size) + " -H -q /dev/stdin -r " + subset_target_fn + \
                      " --coords " + aln_fn
            print(aln_cmd)
            call(aln_cmd, shell=True)

            if not args.keep:
                print("alignment completed; removing the tmp dir.")
                rm_cmd = "rm -rf " + args.tmp_dir
                print(rm_cmd)
                call(rm_cmd, shell=True)

            print("loading jf_aligner's output")
            aln_fh = open(aln_fn, 'r')
            reads_tbl_jf = dict()
            for line in tqdm(aln_fh):
                if line[0] == ">":
                    qname = line.strip().split(" ")[1]
                else:
                    rec = line.strip()
                    if qname not in reads_tbl_jf:
                        reads_tbl_jf[qname] = [rec]
                    else:
                        reads_tbl_jf[qname].append(rec)
            aln_fh.close()

            reads_lst_jf = list(reads_tbl_jf.keys())
            reads_num_jf = len(reads_lst_jf)

            print("building jf_aligner's compatibility matrix")
            cover_tbl, cmpt_tbl_jf = build_tbls_paf(reads_lst_jf, reads_tbl_jf, tx_tbl)

            print("writing out jf_aligner's compatibility matrix")
            cmpt_jf_fn = os.path.join(args.out_dir, "cmpt_tbl_jf.tsv")
            cmpt_jf_fh = open(cmpt_jf_fn, 'w')
            cover_fn = os.path.join(args.out_dir, "cover_tbl.tsv")
            cover_fh = open(cover_fn, 'w')

            for i in tqdm(range(reads_num_jf)):
                qname = reads_lst_jf[i]
                cmpt_jf_fh.write(qname)
                cover_fh.write(qname)
                t_idx_set = cmpt_tbl_jf[i]
                for t_idx in t_idx_set:
                    tname = tx_lst[t_idx]
                    cmpt_jf_fh.write("\t" + tname)
                    cover = cover_tbl[i][t_idx]
                    cover_fh.write("\t(" + tname + ", " + str(cover) + ")")
                cmpt_jf_fh.write("\n")
                cover_fh.write("\n")
            cmpt_jf_fh.close()
            cover_fh.close()

            print("building the score matrix")
            score_tbl = build_score_tbl(cmpt_tbl, reads_lst, cmpt_tbl_jf, reads_lst_jf, ms_tbl, cover_tbl,
                                        args.score_ratio)
            print("writing out the score matrix")
            scores_fn = os.path.join(args.out_dir, "scores.tsv")
            scores_fh = open(scores_fn, 'w')
            for i in tqdm(range(len(reads_lst))):
                qname = reads_lst[i]
                scores_fh.write(qname)
                ti_set = cmpt_tbl[i]
                for ti in ti_set:
                    score = score_tbl[i][ti]
                    tname = tx_lst[ti]
                    scores_fh.write("\t(" + tname + ", " + str(score) + ")")
                scores_fh.write("\n")
            scores_fh.close()


if __name__ == "__main__":
    main()
