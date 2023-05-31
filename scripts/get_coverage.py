# take in the reads file and calculate the gt coverage

import sys
import pyfastx

tx_d = dict()


def load_annot(gtf_fn):
    gtf_fh = open(gtf_fn, 'r')
    global tx_d
    for line in gtf_fh:
        if line[0] == '#':
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            infos = fields[8].split(";")
            for info in infos:
                info = info.strip()
                kv_pair = info.split(" ")
                key = kv_pair[0]
                val = kv_pair[1].replace('"', '')
                if key == "transcript_id":
                    tx_id = val.split(".")[0]
                    if tx_id in tx_d.keys():
                        print("Warning! duplicate transcript_id: " + tx_id)
                    tx_d[tx_id] = 0
                    break
    gtf_fh.close()


def trace_reads(read_fn):
    global tx_d
    fq = pyfastx.Fastx(read_fn)
    for name, seq, _ in fq:
        splits = name.split("_")
        if splits[0] == "unassigned":
            origin = splits[0] + "_" + splits[1] + "_" + splits[2]
        else:
            origin = splits[0] + "_" + splits[1]
        tx_d[origin] += 1


def write_cov(out_fn):
    out_fh = open(out_fn, 'w')
    out_fh.write("transcript_id\tcoverage\n")
    global tx_d
    for tx in tx_d:
        cov = tx_d[tx]
        out_fh.write(tx + "\t" + str(cov) + "\n")


def main(gtf_fn, read_fn, out_fn):
    load_annot(gtf_fn)
    trace_reads(read_fn)
    write_cov(out_fn)
    return None


# 1st arg = reads, 2nd arg = annotation (GTF or GFF)
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
