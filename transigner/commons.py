# commmonly shared imports
from collections import namedtuple
import pysam
from tqdm import tqdm
import pyfastx
import os
import pickle
import sys
import json
from datetime import datetime
import math
import numpy as np
from scipy.special import gammaln
from subprocess import call

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

def split_ln(ln, sep):
    clean_ln = ln.strip()
    fields = clean_ln.split(sep)
    clean_fields = [it.strip() for it in fields]
    return clean_fields

def to_tsv(map, fn):
    lns = []
    for name in map:
        lns.append(f'{name}\t{map[name]}')
    with open(fn, 'w') as f:
        f.write('\n'.join(lns) + '\n')

# TODO: substitute score with ascore
acore = namedtuple('acore', ['rst', 'ren', 'ascore', 'pscore'])
tseg = namedtuple('tseg', ['parent', 'slen', 'sst', 'sen', 'scov', 'si'])