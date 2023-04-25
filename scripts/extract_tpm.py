import sys

tpm_d = dict()


def main(gtf_fn):
    global tpm_d
    fh = open(gtf_fn, 'r')
    cnt = 0
    for line in fh:
        if line[0] == "#":
            continue
        if line.split("\t")[2] == "transcript":
            cols = line.split("\t")[8].split(";")
            rid = "None"
            for col in cols:
                key = col.strip().split(" ")[0].strip()
                val = col.strip().split(" ")[1].strip().replace('"', '')
                if key == "reference_id":
                    cnt += 1
                    rid = val
                elif key == "TPM":
                    tpm = val
                    break
            # None would compile everything without a reference_id
            tpm_d[rid] = tpm
    out_fh = open("./transcripts_tpm.tsv", 'w')
    out_fh.write("transcript_id\ttpm\n")
    print(cnt)
    for rid in tpm_d.keys():
        tpm = tpm_d[rid]
        if rid != "None":
            out_fh.write(rid + "\t" + tpm + "\n")
    out_fh.close()


if __name__ == "__main__":
    main(sys.argv[1])