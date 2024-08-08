RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

def split_ln(ln, sep):
    clean_ln = ln.strip()
    fields = clean_ln.split(sep)
    clean_fields = [it.strip() for it in fields]
    return clean_fields