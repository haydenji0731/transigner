from datetime import datetime
import sys
import os
import json
from enum import Enum
from subprocess import call
import argparse

RED = '\033[31m'
GREEN = '\033[32m'
RESET = '\033[0m'

class Mtype(Enum):
    PROG = (GREEN, "PROGRESS")
    ERR = (RED, "ERROR")
    WARN = (RED, "WARNING")

def store_params(args, fn):
    with open(fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

def tmessage(s, mtype) -> str:
    if mtype not in Mtype:
        raise Exception("Error while printing message")
    return f"{datetime.now()} {mtype.value[0]}{mtype.value[1]}{RESET} {s}"

def check_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)
