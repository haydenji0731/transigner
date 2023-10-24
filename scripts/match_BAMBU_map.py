#!/usr/bin/env python

import sys


def load_tx_tbl(fn):
    fh = open(fn, 'r')
    first = True
    tx_tbl = dict()
    for line in fh:
        if first:
            first = False
            continue
        fields = line.split("\t")
        assert len(fields) == 9
        tid = fields[1].replace('"', '')
        idx = int(fields[0])
        if tid in tx_tbl:
            assert tx_tbl[idx] == tid
        else:
            tx_tbl[idx] = tid
    fh.close()
    return tx_tbl


def load_mapping(fn, tx_tbl):
    fh = open(fn, 'r')
    first = True
    read2tx_tbl = dict()
    for line in fh:
        line = line.strip()
        if first:
            first = False
            continue
        fields = line.split("\t")
        assert len(fields) == 3
        rid = fields[0].replace('"', '')
        eq_match = fields[1].replace('"', '')
        cmpt_match = fields[2].replace('"', '')
        if eq_match == "NULL" and cmpt_match == "NULL":
            read2tx_tbl[rid] = None
        else:
            read2tx_tbl[rid] = list()
            if eq_match != "NULL":
                eq_match_lst = proc_match(eq_match)
                for eq_idx in eq_match_lst:
                    tid = tx_tbl[eq_idx]
                    # 0 means equal match
                    read2tx_tbl[rid].append((tid, 0))
            if cmpt_match != "NULL":
                cmpt_match_lst = proc_match(cmpt_match)
                for cmpt_idx in cmpt_match_lst:
                    tid = tx_tbl[cmpt_idx]
                    # 1 means compatible match
                    read2tx_tbl[rid].append((tid, 1))
    return read2tx_tbl


# process a line containing match info (compatible or equal)
def proc_match(ln):
    match_lst = list()
    if ln[0] == 'c':
        tmp = ln.split(",")
        for i in range(len(tmp)):
            idx = tmp[i]
            if i == 0:
                idx = int(idx.replace('c(', ''))
            elif i == len(tmp) - 1:
                idx = int(idx.replace(')', ''))
            else:
                idx = int(idx.strip())
            match_lst.append(idx)
    else:
        tmp = ln.split(":")
        if len(tmp) > 1:
            st = int(tmp[0])
            en = int(tmp[1])
            for i in range(st, en + 1):
                match_lst.append(i)
        else:
            idx = int(tmp[0])
            match_lst.append(idx)
    return match_lst


def main(tx_fn, map_fn, out_fn):
    tx_tbl = load_tx_tbl(tx_fn)
    read2tx_tbl = load_mapping(map_fn, tx_tbl)
    out_fh = open(out_fn, 'w')
    for rid in read2tx_tbl:
        match_lst = read2tx_tbl[rid]
        out_fh.write(rid)
        if match_lst is None:
            out_fh.write("\tNone\n")
            continue
        for match in match_lst:
            out_fh.write("\t" + match[0])
        out_fh.write("\n")
    out_fh.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])



