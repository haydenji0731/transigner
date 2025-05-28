#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import os
import sys

def load_counts(fn):
    df = pd.read_csv(fn, sep=',', header=None)
    df.columns = ['transcript_id', 'read_count']
    return df

def load_assignments(fn):
    assignments = dict()
    with open(fn, 'r') as fh:
        for ln in fh:
            parts = ln.split("\t")
            for i in range(len(parts)):
                if i == 0:
                    read_id = parts[i]
                    assignments[read_id] = []
                else:
                    temp = parts[i].split(',')
                    assignments[read_id].append((temp[0], float(temp[1])))
    return assignments
                
def calc_rmse(y, y_hat):
    mse = ((y - y_hat) ** 2).mean()
    return np.sqrt(mse)

def get_gt(s):
    temp = s.split('_')
    if temp[0] == "unassigned":
        gt = '_'.join(temp[0:3])
    else:
        gt = '_'.join(temp[0:2])
    return gt

def calc_pr(assignments, n):
    tp = 0
    for read_id in assignments:
        gt_transcript_id = get_gt(read_id)
        for transcript_id, frac in assignments[read_id]:
            if transcript_id.split('.')[0] == gt_transcript_id:
                tp += frac
    recall = tp / n
    precision = tp / len(assignments)
    return precision, recall

def main(args):
    y_df = load_counts(args.y)
    y_hat_df = load_counts(os.path.join(args.y_hat, 'abundances.csv'))
    merged_y = pd.merge(y_df, y_hat_df, on='transcript_id', how='left')
    merged_y = merged_y.fillna(0)

    pearson_r = np.log2(merged_y['read_count_x'] + 1).corr(np.log2(merged_y['read_count_y'] + 1), method='pearson')
    spearman_r = merged_y['read_count_x'].corr(merged_y['read_count_y'], method='spearman')
    rmse = calc_rmse(merged_y['read_count_x'], merged_y['read_count_y'])

    assignments = load_assignments(os.path.join(args.y_hat, 'assignments.tsv'))
    precision, recall = calc_pr(assignments, int(args.n))

    print(f"pearson's\t{pearson_r:.3f}\nspearman's\t{spearman_r:.3f}\nrmse\t{rmse:.3f}")
    print(f"precision\t{precision:.3f}\nrecall\t{recall:.3f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-y')
    parser.add_argument('-y-hat')
    parser.add_argument('-n')
    args = parser.parse_args()
    if not os.path.isdir(args.y_hat):
        print("Error: --y-hat dir not present")
        sys.exit(-1)
    main(args)
    

