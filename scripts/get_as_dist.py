#!/usr/bin/env python

import pysam
import sys
from time import strftime

primary = dict()
secondary = dict()
score_dist = dict()


def load_aln(aln_fn):
    unmapped = 0
    unique = 0
    global primary
    global secondary
    qname_prev = ""
    secondary_l = list()
    with pysam.AlignmentFile(aln_fn, 'rb') as fh:
        for brec in fh:
            if brec.is_unmapped:
                unmapped += 1
            else:
                qname_curr = brec.query_name
                if qname_prev != qname_curr:
                    # if not 1st read
                    if qname_prev != "":
                        secondary[qname_prev] = secondary_l
                        if len(secondary_l) == 0:
                            unique += 1
                        secondary_l = list()
                    qname_prev = qname_curr
                # chain_score = int(brec.get_tag("s1"))
                # num of mismatches + gaps
                # nm_score = int(brec.get_tag("NM"))
                # max scoring segment
                ms_score = int(brec.get_tag("ms"))
                if not brec.is_secondary and not brec.is_supplementary:
                    primary[qname_curr] = ms_score
                elif brec.is_secondary:
                    secondary_l.append(ms_score)
    if len(secondary_l) == 0:
        unique += 1
    secondary[qname_prev] = secondary_l
    # for sanity check (what tag is used to determine primary alignment?)
    for read in primary.keys():
        pri_score = primary[read]
        for sa in secondary[read]:
            if pri_score < sa:
                print("anomaly detected using ms score")
    return unmapped, unique


def get_score_dist(out_fh, out_raw_fh):
    global primary
    global secondary
    pri_zero = 0
    sec_eq = 0
    sec_gt = 0
    sec_lt = 0
    unique = 0
    for qname in primary.keys():
        status = None
        pri_score = primary[qname]
        sec_l = secondary[qname]
        sec_l.sort()
        # focus on multi-mapped reads for now
        sec_frac_l = list()
        if pri_score == 0:
            print("primary alignment score equals 0; skipping read: " + qname)
            pri_zero += 1
            continue
        if len(sec_l) > 0:
            for sec in sec_l:
                s2p_ratio = sec / pri_score
                sec_frac = (sec - pri_score) / abs(pri_score)
                if status is None:
                    if sec_frac == 0:
                        status = "eq"
                    elif sec_frac > 0:
                        status = "gt"
                    else:
                        status = "lt"
                sec_frac_l.append(sec_frac)
            sec_min = min(sec_frac_l)
            sec_max = max(sec_frac_l)
            sec_mean = sum(sec_frac_l) / len(sec_frac_l)
            out_fh.write(qname + "\t" + str(pri_score) + "\t" + str(sec_min) + "\t" + str(sec_max)
                         + "\t" + str(sec_mean) + "\t(")
            for frac in sec_frac_l:
                if frac == 121.2:
                    print("max: " + qname)
                if frac == -14.25:
                    print("min: " + qname)
                out_fh.write(str(frac) + ",")
                out_raw_fh.write(str(frac) + "\n")
            out_fh.write(")\n")
        if status == "eq":
            sec_eq += 1
        elif status == "gt":
            sec_gt += 1
        elif status == "lt":
            sec_lt += 1
        else:
            unique += 1
    print("# of reads with secondary scores equal to primary: " + str(sec_eq))
    print("# of reads with secondary scores greater than primary: " + str(sec_gt))
    print("# of reads with secondary scores less than primary: " + str(sec_lt))
    print("# of reads with 0 primary score: " + str(pri_zero))
    print("# of unique: " + str(unique))


def main(aln_fn, out_fn, out_raw_fn):
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading alignment file")
    unmapped, unique = load_aln(aln_fn)
    print("# of unmapped reads: " + str(unmapped))
    print("# of uniquely-mapped reads: " + str(unique))
    out_fh = open(out_fn, 'w')
    out_fh.write("read_id\tpri\tsec_min\tsec_max\tsec_mean\tsec_raw\n")
    out_raw_fh = open(out_raw_fn, 'w')
    out_raw_fh.write("secondary_frac\n")
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Obtaining score distribution")
    get_score_dist(out_fh, out_raw_fh)
    out_fh.close()
    out_raw_fh.close()
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished writing distributions to the output file")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
