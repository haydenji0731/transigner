#!/usr/bin/env python

import sys
import re


def main(fn, out_fn):
    fh = open(fn, 'r')
    out_fh = open(out_fn, 'w')
    pat1 = "c\("
    pat2 = "\)"
    truncated = False

    for line in fh:
        if truncated:
            res = re.search(pat2, line)
            if res is None:
                tmp += line.strip()
            else:
                tmp += line
                out_fh.write(tmp)
                truncated = False
            continue
        res = re.search(pat1, line)
        if res is None:
            out_fh.write(line)
        else:
            res = re.search(pat2, line)
            if res is None:
                truncated = True
                tmp = line.strip()
    fh.close()
    out_fh.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
