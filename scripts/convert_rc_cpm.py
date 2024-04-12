#!/usr/bin/env python

import sys
import math

def main(in_fn, ls, out_fn):
	library_size = int(ls) / 4
	in_fh = open(in_fn, 'r')
	out_fh = open(out_fn, 'w')
	for ln in in_fh:
		clean_ln = ln.strip()
		fields = clean_ln.split("\t")
		rc = float(fields[1])
		cpm = rc / library_size * 1000000
		log_cpm = math.log2(cpm + 1)
                # converted this to output cpm instead
		out_fh.write(fields[0] + "\t" + str(cpm) + "\n")
	in_fh.close()
	out_fh.close()


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])
