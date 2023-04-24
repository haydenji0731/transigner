#!/usr/bin/env python

import sys

gtf_d = dict()
tpm_d = dict()


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
                rid = None
                if key == "transcript_id":
                    tid = val
                elif key == "reference_id":
                    rid = val.split(".")[0]
                elif key == "TPM":
                    tpm = float(val)
                    break
            gtf_d[tid] = (rid, tpm)


def load_true_tpm(fname):
    fh = open(fname, 'r')
    first = True
    global tpm_d
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        tpm = float(line.split("\t")[2])
        tpm_d[tid] = tpm


def calc_loss(fname):
    fh = open(fname, 'r')
    first = True
    st_loss = 0
    loss = 0
    global gtf_d
    global tpm_d
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        tpm = float(line.split("\t")[3])
        rid = gtf_d[tid][0]
        st_tpm = gtf_d[tid][1]
        if rid is not None:
            true_tpm = tpm_d[rid]
            tpm_d[rid] = 0
            loss += abs(true_tpm - tpm)
            st_loss += abs(true_tpm - st_tpm)
        else:
            loss += tpm
            st_loss += st_tpm
        # remaining = not assembled but positive TPM
    remaining = sum(tpm_d.values())
    loss += remaining
    st_loss += remaining

    # total number of reference transcripts
    tot = len(tpm_d.keys())
    loss /= tot
    st_loss /= tot
    return loss, st_loss


def main(gtf_fn, label_fn, abundance_fn):
    load_gtf(gtf_fn)
    load_true_tpm(label_fn)
    loss, st_loss = calc_loss(abundance_fn)
    print("StringTie2 mean absolute loss: " + str(st_loss))
    print("TranSigner mean absolute loss: " + str(loss))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])