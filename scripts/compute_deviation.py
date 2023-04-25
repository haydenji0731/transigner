#!/usr/bin/env python

import sys

gtf_d = dict()


def load_gtf(fname):
    global gtf_d
    fh = open(fname, 'r')
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
                elif key == "TPM":
                    tpm = float(val)
                    break
            gtf_d[tid] = tpm


def calc_deviation(fname):
    global gtf_d
    fh = open(fname, 'r')
    first = True
    sigma = 0
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        tpm = float(line.split("\t")[3])
        gtf_tpm = gtf_d[tid]
        sigma += abs(gtf_tpm - tpm)
    tot = len(gtf_d.keys())
    sigma /= tot
    return sigma


def main(gtf_fn, abundance_fn):
    load_gtf(gtf_fn)
    sigma = calc_deviation(abundance_fn)
    print("mean deviation: " + str(sigma))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])