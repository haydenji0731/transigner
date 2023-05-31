import sys


def main(gtf_fn, out_fn):
    gtf_fh = open(gtf_fn, 'r')
    out_fh = open(out_fn, 'w')
    for line in gtf_fh:
        if line[0] == '#':
            continue
        fields = line.split("\t")
        infos = fields[8].split(";")
        for info in infos:
            clean_info = info.strip()
            kv_pair = clean_info.split(" ")
            key = kv_pair[0]
            val = kv_pair[1].replace('"', '')
            if key == "transcript_id":
                tx_id = val.replace('.', '_')
                break
        for i in range(len(fields)):
            if i == 8:
                for info in infos:
                    if info == ' \n':
                        out_fh.write(info)
                        continue
                    clean_info = info.strip()
                    kv_pair = clean_info.split(" ")
                    key = kv_pair[0]
                    if key == "transcript_id":
                        out_fh.write(" transcript_id " + '"' + tx_id + '";')
                    else:
                        out_fh.write(info + ";")
            else:
                out_fh.write(fields[i] + "\t")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])