from transigner.commons import split_ln
import sys

class Line():
    def __init__(self, ln, format):
        if ln is None:
            self.init_empty()
        else:
            if format == "gtf":
                kv_sep = ' '
                quoted = True
            elif format == "gff":
                kv_sep = '='
                quoted = False
            else:
                print("Invalid format info")
                sys.exit(-1)
            
            fields = split_ln(ln, sep="\t")

            if len(fields) != 9:
                print("A line must have exactly 9 columns.")
                sys.exit(-1)

            self.ctg = fields[0]
            self.src = fields[1]
            self.feature = fields[2]
            self.start = int(fields[3])
            self.end = int(fields[4])
            self.score = float(fields[5]) if fields[5] != '.' else None
            self.strand = fields[6]
            self.frame = fields[7] if fields[7] != '.' else None
            
            # split attributes into a dictionary
            tmp_arr = split_ln(fields[8], sep=';')
            self.attributes = dict()
            for att in tmp_arr:
                att_kv_pair = att.split(kv_sep)
                if len(att_kv_pair) < 2:
                    continue
                k = att_kv_pair[0]
                v = att_kv_pair[1]
                if quoted:
                    v = v.replace('"', '')
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
    
    def get_feature_id(self, att_pattern):
        # print(self.attributes)
        feature_id = None
        for att in self.attributes:
            if att == att_pattern:
                feature_id = self.attributes[att]
        return feature_id
