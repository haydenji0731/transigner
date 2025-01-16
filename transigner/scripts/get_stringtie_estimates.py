#!/usr/bin/env python

import argparse
import os

class Line():
    def __init__(self, ln):
        if ln is None:
            self.init_empty()
        fields = [x.strip() for x in ln.strip().split('\t')]
        if len(fields) != 9: raise Exception("An input line must have exactly 9 columns"); sys.exit(-1)
        self.ctg = fields[0]
        self.src = fields[1]
        self.feature = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.score = float(fields[5]) if fields[5] != '.' else None
        self.strand = fields[6]
        self.frame = fields[7] if fields[7] != '.' else None
        self.attributes = dict()
        for att in [x.strip() for x in fields[8].split(';')]:
            att_kv_pair = att.split(" ")
            if len(att_kv_pair) < 2: continue
            k = att_kv_pair[0]
            v = att_kv_pair[1].replace('"', '')
            self.attributes[k] = v

    def init_empty(self):
        self.ctg = None
        self.src = None
        self.feature = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.frame = None
        self.attributes = None
        
def load_basic(in_fn):
    estimates = dict()
    with open(in_fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_obj = Line(ln)
            if ln_obj.feature == "transcript":
                tid = ln_obj.attributes['transcript_id']
                cov = ln_obj.attributes['cov']
                fpkm = ln_obj.attributes['FPKM']
                tpm = ln_obj.attributes['TPM']
                estimates[tid] = (cov, fpkm, tpm)
    return estimates

def load_full(in_fn):
    tx_estimates = dict()
    exon_estimates = dict()
    with open(in_fn, 'r') as fh:
        for ln in fh:
            if ln[0] == '#': continue
            ln_obj = Line(ln)
            tid = ln_obj.attributes['transcript_id']
            if ln_obj.feature == "transcript":
                cov = float(ln_obj.attributes['cov'])
                fpkm = float(ln_obj.attributes['FPKM'])
                tpm = float(ln_obj.attributes['TPM'])
                tx_estimates[tid] = (cov, fpkm, tpm, ln_obj.strand)
            elif ln_obj.feature == "exon":
                cov = float(ln_obj.attributes['cov'])
                if tid in exon_estimates:
                    # TODO: maybe lose the start and end coords?
                    exon_estimates[tid].append((ln_obj.start, ln_obj.end, cov))
                else:
                    exon_estimates[tid] = [(ln_obj.start, ln_obj.end, cov)]
    return tx_estimates, exon_estimates

def print_transcript_estimates(e, fn):
    with open(fn, 'w') as fh:
        fh.write(f'transcript_id,cov,fpkm,tpm\n')
        for x in e: temp = e[x]; fh.write(f'{x},{temp[0]},{temp[1]},{temp[2]}\n')

def print_exon_estimates(e, fn):
    with open(fn, 'w') as fh:
        fh.write(f'transcript_id,exon_id,start,end,cov\n')
        for x in e:
            for i, ex in enumerate(e[x]): fh.write(f'{x},{x}-exon-{i+1},{ex[0]},{ex[1]},{ex[2]}\n')

def main() -> None:
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--in-file', type=str, help="", required=True)
    parser.add_argument('-o', '--out-dir', type=str, help="", required=True)
    parser.add_argument('--include-exons', default=False, help="", action='store_true')
    args = parser.parse_args()
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    if args.include_exons:
        tx_e, exon_e = load_full(args.in_file)
        out_fn = os.path.join(args.out_dir, 'transcript_estimates.csv')
        print_transcript_estimates(tx_e, out_fn)
        out_fn = os.path.join(args.out_dir, 'exon_estimates.csv')
        print_exon_estimates(exon_e, out_fn)
    else:
        e = load_basic(args.in_file)
        out_fn = os.path.join(args.out_dir, 'transcript_estimates.csv')
        print_transcript_estimates(e, out_fn)

if __name__ == "__main__":
    main()