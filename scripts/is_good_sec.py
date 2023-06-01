#!/usr/bin/env python

import sys
import pysam
from time import strftime

ref2tx = dict()
ambig = set()


def get_gt(gtf_fn):
    gtf_fh = open(gtf_fn, 'r')
    global ref2tx
    global ambig
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
                elif key == "reference_id":
                    ref_id = val.split('.')[0]
                    break
            if not no_ref:
                if ref_id in ref2tx.keys():
                    # TODO: overcome the NanoSim menace
                    print("Warning: > 1 transcript detected for the same reference transcript")
                    ambig.add(ref_id)
                else:
                    ref2tx[ref_id] = tx_id


def match_gt2aln(aln_fn, out_fn, absent_fn):
    global ref2tx
    global ambig
    absent = set()
    qname_prev = ""
    unassembled_cnt = 0
    sec_match_gt_wls = 0
    out_fh = open(out_fn, 'w')
    out_fh.write("read_id\tprimary_score\tsecondary_score\n")
    with pysam.AlignmentFile(aln_fn, 'rb') as fh:
        for brec in fh:
            if not brec.is_unmapped:
                qname = brec.query_name
                if qname != qname_prev:
                    if qname_prev != "":
                        if is_sec_gt:
                            if pri_score is None:
                                raise Exception("ERROR: no primary score detected for this read")
                            if sec_score < pri_score:
                                sec_match_gt_wls += 1
                                out_fh.write(qname + "\t" + str(pri_score) + "\t" + str(sec_score) + "\n")
                    pri_score = None
                    sec_score = None
                    is_sec_gt = False
                    qname_prev = qname
                    unassembled = False
                if unassembled:
                    continue

                if not brec.is_secondary and not brec.is_supplementary:
                    pri_score = int(brec.get_tag("ms"))
                # if not is_sec_gt:
                qname_pieces = qname.split("_")
                if qname_pieces[0] == "unassigned":
                    ref_id = qname_pieces[0] + "_" + qname_pieces[1] + "_" + qname_pieces[2]
                else:
                    ref_id = qname_pieces[0] + "_" + qname_pieces[1]
                if ref_id in ambig:
                    print("Warning: an ambiguous reference id encountered!")
                if ref_id in absent:
                    unassembled = True
                elif ref_id not in ref2tx.keys():
                    # print("Warning: this transcript wasn't assembled - " + ref_id)
                    unassembled_cnt += 1
                    unassembled = True
                    absent.add(ref_id)
                else:
                    gt = ref2tx[ref_id]
                    tid = brec.reference_name
                    if gt == tid:
                        if brec.is_secondary:
                            is_sec_gt = True
                            sec_score = int(brec.get_tag("ms"))

    if is_sec_gt:
        if pri_score is None:
            raise Exception("ERROR: no primary score detected for this read")
        if sec_score < pri_score:
            sec_match_gt_wls += 1
            out_fh.write(qname + "\t" + str(pri_score) + "\t" + str(sec_score) + "\n")
    absent_fh = open(absent_fn, 'w')
    absent_fh.write("transcript_id\n")
    for tx in absent:
        absent_fh.write(tx + "\n")
    return unassembled_cnt, sec_match_gt_wls


def main(aln_fn, gtf_fn, out_fn, absent_fn):
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Matching ground truth to assembled transcriptome")
    get_gt(gtf_fn)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Checking for good secondary alignments")
    unassembled_cnt, good_sec = match_gt2aln(aln_fn, out_fn, absent_fn)
    print("# of unassembled ground truth transcripts: " + str(unassembled_cnt))
    print("# of good secondary alignments: " + str(good_sec))


if __name__ == "__main__":
    # input BAM, input gtf
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
