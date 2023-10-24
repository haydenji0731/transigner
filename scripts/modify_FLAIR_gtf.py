#!/usr/bin/env python

import sys


def mod_gtf(fn, out_fn):
    fh = open(fn, 'r')
    out_fh = open(out_fn, 'w')
    for line in fh:
        fields = line.split('\t')
        for i in range(len(fields)):
            if i < 8:
                out_fh.write(fields[i] + "\t")
            else:
                infos = fields[i].split(";")
                infos = infos[:-1]
                for j in range(len(infos)):
                    info = infos[j]
                    clean_info = info.strip()
                    kv_pair = clean_info.split(" ")
                    k = kv_pair[0]
                    v = kv_pair[1].replace('"', '')
                    if k == "gene_id":
                        gid = v
                        out_fh.write('gene_id "' + gid + '"; ')
                    elif k == "transcript_id":
                        tid = v
                        mod_tid = tid + "_" + gid
                        if j == len(infos) - 1:
                            out_fh.write('transcript_id "' + mod_tid + '";\n')
                        else:
                            out_fh.write('transcript_id "' + mod_tid + '"; ')
                    elif j == len(infos) - 1:
                        out_fh.write(k + ' "' + v + '";\n')
                    else:
                        out_fh.write(k + ' "' + v + '"; ')
    fh.close()
    out_fh.close()


def main(gtf_fn, out_fn):
    mod_gtf(gtf_fn, out_fn)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

