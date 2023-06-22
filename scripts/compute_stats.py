#!/usr/bin/env python

import sys
import pysam
from time import strftime

ref2tx = dict()
aln_d = dict()
dup = set()


def get_gt(gtf_fn, out_prefix):
    gtf_fh = open(gtf_fn, 'r')
    out_fn = out_prefix + "_gt.tsv"
    out_fh = open(out_fn, 'w')
    global ref2tx
    global dup
    for line in gtf_fh:
        if line[0] == '#':
            continue
        feature = line.split("\t")[2]
        if feature == "transcript":
            no_ref = False
            infos = line.split("\t")[8].split(";")
            for info in infos:
                if info == '\n':
                    no_ref = True
                    break
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tx_id = val
                elif key == "cmp_ref":
                    # TODO: fix this since the reads are also fixed
                    ref_id = val.split('.')[0]
                    # ref_id = val
                elif key == "class_code":
                    code = val
                    break
            if not no_ref:
                if code == '=':
                    if ref_id in ref2tx.keys():
                        print("duplicate detected")
                        print(ref_id + "\t" + tx_id + "\n")
                        dup.add(ref_id)
                        ref2tx[ref_id].add(tx_id)
                    else:
                        ref2tx[ref_id] = {tx_id}
            else:
                print("Warning: no reference found for transcript id: " + tx_id)
    for ref in ref2tx.keys():
        out_fh.write(ref + "\t")
        tx_l = list(ref2tx[ref])
        for i in range(len(tx_l)):
            if i != len(tx_l) - 1:
                out_fh.write(tx_l[i] + "\t")
            else:
                out_fh.write(tx_l[i] + "\n")
    print("total # of reference transcripts with multiple matches: " + str(len(dup)))

    # OLD code for using cmp_ref field of guided StringTie assembly --> previous approach (potentially incorrect)
    # if not no_ref:
    #     if ref_id in ref2tx.keys():
    #         print("Warning: > 1 transcript detected for the same reference transcript")
    #         ambig.add(ref_id)
    #     else:
    #         ref2tx[ref_id] = tx_id


# def find_gt_from_name(qname):
#     qname_pieces = qname.split("_")
#     gt_ref_id = qname_pieces[0] + "_" + qname_pieces[1] + "_" + qname_pieces[2]
#     return gt_ref_id

def find_gt_from_name(qname):
    qname_pieces = qname.split("_")
    if qname_pieces[0] == "unassigned":
        gt_ref_id = qname_pieces[0] + "_" + qname_pieces[1] + "_" + qname_pieces[2]
    else:
        gt_ref_id = qname_pieces[0] + "_" + qname_pieces[1]
    return gt_ref_id


def load_aln(aln_fn):
    global aln_d
    unmapped = 0
    with pysam.AlignmentFile(aln_fn, 'rb') as fh:
        for brec in fh:
            if brec.is_unmapped:
                unmapped += 1
            else:
                qname = brec.query_name
                if qname not in aln_d.keys():
                    aln_d[qname] = [brec]
                else:
                    aln_d[qname].append(brec)
    mapped = len(aln_d.keys())
    return unmapped, mapped


def compute_stats(aln_fn, out_prefix, unmapped, mapped):
    # only contains MAPPED reads
    global aln_d
    missing_gt = 0
    unassembled = set()
    multi = 0
    unique = 0
    w_supl = 0
    good_sec_wls = 0
    good_sec = 0
    good_supl = 0
    scores_fn = out_prefix + "_scores.tsv"
    stats_fn = out_prefix + "_stats.tsv"
    missing_fn = out_prefix + "_missing.txt"
    scores_fh = open(scores_fn, 'w')
    stats_fh = open(stats_fn, 'w')
    missing_fh = open(missing_fn, 'w')
    scores_fh.write("read_id\ttarget_id\tscore\ttype\tnote\n")
    stats_fh.write("Stats for input file: " + aln_fn + "\n\n")
    stats_fh.write("================================================\n")
    for qname in aln_d.keys():
        gt_ref_id = find_gt_from_name(qname)
        has_sec = False
        has_supl = False
        is_sec_lt = False
        is_sec_good = False
        is_supl_good = False
        sec_l = list()
        supl_l = list()
        pri_score = None
        pri_tid = None
        if gt_ref_id not in ref2tx.keys():
            missing_gt += 1
            unassembled.add(gt_ref_id)
            # TODO: save unassembled gt ids and output
            continue
        else:
            brecs = aln_d[qname]
            for brec in brecs:
                if brec.is_secondary:
                    has_sec = True
                    target_id = brec.reference_name
                    ref_ids = ref2tx[gt_ref_id]
                    if target_id in ref_ids:
                        is_sec_good = True
                        sec_score = int(brec.get_tag("ms"))
                        sec_l.append((target_id, sec_score))
                elif brec.is_supplementary:
                    has_supl = True
                    target_id = brec.reference_name
                    ref_ids = ref2tx[gt_ref_id]
                    if target_id in ref_ids:
                        is_supl_good = True
                        supl_score = int(brec.get_tag("ms"))
                        supl_l.append((target_id, supl_score))
                elif not brec.is_secondary and not brec.is_supplementary:
                    pri_score = int(brec.get_tag("ms"))
                    pri_tid = brec.reference_name
                else:
                    print("WARNING: this is an unhandled edge case")
            if has_sec:
                multi += 1
            if has_supl:
                w_supl += 1
            if not has_sec or not has_supl:
                unique += 1
            # output
            if pri_score is None and pri_tid is None:
                print("WARNING: no primary detected for read: " + qname)
                continue
            else:
                scores_fh.write(qname + "\t" + pri_tid + "\t" + str(pri_score) +
                                "\tPR\t.\n")
            if is_sec_good:
                good_sec += 1
                for sec in sec_l:
                    if sec[1] < pri_score:
                        is_sec_lt = True
                        scores_fh.write(qname + "\t" + sec[0] + "\t" + str(sec[1]) +
                                        "\tSC\twls:1\n")
                    else:
                        scores_fh.write(qname + "\t" + sec[0] + "\t" + str(sec[1]) +
                                        "\tSC\twls:0\n")
            if is_sec_lt:
                good_sec_wls += 1
            if is_supl_good:
                good_supl += 1
                for supl in supl_l:
                    scores_fh.write(qname + "\t" + supl[0] + "\t" + str(supl[1]) +
                                    "\tSP\t.\n")
    tot = mapped + unmapped
    stats_fh.write("%d / %d reads mapped\n" % (mapped, tot))
    stats_fh.write("%d / %d reads unmapped\n" % (unmapped, tot))
    stats_fh.write("%d / %d reads whose ground truth wasn't assembled\n" % (missing_gt, tot))
    stats_fh.write("%d / %d reads multi-mapped\n" % (multi, tot))
    stats_fh.write("%d / %d reads uniquely mapped\n" % (unique, tot))
    stats_fh.write("%d / %d reads are chimeric\n" % (w_supl, tot))
    stats_fh.write("%d reads with good secondary alignments\n" % good_sec)
    stats_fh.write("%d / %d reads with good secondary alignments has < score than primary\n" %
                   (good_sec_wls, good_sec))
    stats_fh.write("%d reads with good supplementary alignments\n" % good_supl)

    for rid in unassembled:
        missing_fh.write(rid + "\n")

    scores_fh.close()
    stats_fh.close()
    missing_fh.close()


def main(aln_fn, gtf_fn, out_prefix):
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Matching ground truth to assembled transcriptome")
    get_gt(gtf_fn, out_prefix)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading input alignment file")
    unmapped, mapped = load_aln(aln_fn)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Calculating stats around good secondary and supplementary alignments")
    compute_stats(aln_fn, out_prefix, unmapped, mapped)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Processing completed")


if __name__ == "__main__":
    # input BAM, input gtf, output prefix
    main(sys.argv[1], sys.argv[2], sys.argv[3])
