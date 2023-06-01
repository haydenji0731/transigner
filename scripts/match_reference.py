import sys

gtf_d = dict()


def load_gtf(fname):
    global gtf_d
    fh = open(fname, 'r')
    for line in fh:
        if line[0] == "#":
            continue
        if line.split("\t")[2] == "transcript":
            cols = line.split("\t")[8].split(";")
            for col in cols:
                key = col.strip().split(" ")[0].strip()
                val = col.strip().split(" ")[1].strip().replace('"', '')
                if key == "transcript_id":
                    tid = val
                elif key == "reference_id":
                    rid = val
                    break
                elif key == "cov":
                    rid = "None"
                    break
            gtf_d[tid] = rid


def main(gtf_fn, abundance_fn):
    load_gtf(gtf_fn)
    fh = open(abundance_fn, 'r')
    out_fh = open("./abundance_modified.tsv", 'w')
    out_fh.write("transcript_id\ttpm\n")
    first = True
    global gtf_d
    for line in fh:
        if first:
            first = False
            continue
        tid = line.split("\t")[0]
        tpm = line.split("\t")[3]
        rid = gtf_d[tid]
        if rid != "None":
            out_fh.write(rid + "\t" + tpm)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
