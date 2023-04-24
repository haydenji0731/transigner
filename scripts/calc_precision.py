#!/usr/bin/env python

import sys

gtf_d = dict()
# tn = unaligned
# fn = aligned but unmapped


def load_gtf(fname):
    global gtf_d
    fh = open(fname, 'r')
    novel = 0
    for line in fh:
        if line[0] == "#":
            continue
        if line.split("\t")[2] == "transcript":
            infos = line.split("\t")[8].split(";")
            for info in infos:
                key = info.strip().split(" ")[0].strip()
                val = info.strip().split(" ")[1].strip().replace('"', '')
                if key == "transcript_id":
                    tid = val
                elif key == "reference_id":
                    rid = val.split(".")[0]
                    break
                elif key == "cov":
                    rid = None
                    novel += 1
                    break
        gtf_d[tid] = rid
    return novel


def load_assignment(fname):
    fh = open(fname, 'r')
    first = True
    tp = 0
    fp = 0
    tp_conf = 0
    fp_conf = 0
    global gtf_d
    for line in fh:
        if first:
            first = False
            continue
        tmp = line.split("\t")[0].split("_")
        pred = tmp[0] + tmp[1]
        tid = line.split("\t")[1]
        label = gtf_d[tid]
        conf = float(line.split("\t")[2])
        if label is None or pred != label:
            fp += 1
            fp_conf += conf
        elif pred == label:
            tp += 1
            tp_conf += conf
    # take average
    fp_conf /= fp
    tp_conf /= tp
    return tp, tp_conf, fp, fp_conf


def main(gtf_fn, assignment_fn):
    novel = load_gtf(gtf_fn)
    print("# of novel transcripts: " + str(novel))
    tp, tp_conf, fp, fp_conf = load_assignment(assignment_fn)
    print("# of True Positives: " + str(tp))
    print("average confidence for True Positives: " + str(tp_conf))
    print("# of False Positives: " + str(fp))
    print("average confidence for False Positives: " + str(fp_conf))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])