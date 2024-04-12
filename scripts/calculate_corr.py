#!/usr/bin/env python

import sys
import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
import math


def compute_ci(r, n, ci=0.95):
    z_score = st.norm.ppf((1 + ci) / 2)
    stderr = 1 / math.sqrt(n - 3)
    delt = z_score * stderr
    lo = math.tanh(math.atanh(r) - delt)
    hi = math.tanh(math.atanh(r) + delt)
    return lo, hi


def load_cnt(fn):
    cnt_tbl = dict()
    fh = open(fn, 'r')
    for line in fh:
        clean_ln = line.strip()
        tmp = clean_ln.split("\t")
        tid = tmp[0]
        cnt = float(tmp[1])
        cnt_tbl[tid] = cnt
    return cnt_tbl


def extract_cnt(gt_cnt_tbl, est_cnt_tbl, is_linear):
    gt_cnt_lst = list()
    est_cnt_lst = list()
    for tid in gt_cnt_tbl:
        gt_cnt = gt_cnt_tbl[tid]
        if is_linear == 0:
            gt_cnt = np.log2(gt_cnt + 1)
        gt_cnt_lst.append(gt_cnt)
        est_cnt = est_cnt_tbl[tid]
        if is_linear == 0:
            est_cnt = np.log2(est_cnt + 1)
        est_cnt_lst.append(est_cnt)
    assert len(gt_cnt_lst) == len(est_cnt_lst)
    return gt_cnt_lst, est_cnt_lst


def main(gt_fn, est_fn, save_fn, is_linear, is_svg):
    gt_cnt_tbl = load_cnt(gt_fn)
    tmp_tbl = load_cnt(est_fn)
    est_cnt_tbl = dict()
    for tid in tmp_tbl:
        if tid in gt_cnt_tbl:
            est_cnt_tbl[tid] = tmp_tbl[tid]
    gt_cnt_lst, est_cnt_lst = extract_cnt(gt_cnt_tbl, est_cnt_tbl, is_linear)

    # sanity check
    assert len(gt_cnt_lst) == len(est_cnt_lst)

    if is_linear == 0:
        spearman_res = st.spearmanr(gt_cnt_lst, est_cnt_lst)
        spearman_r = spearman_res.statistic
        spearman_pval = spearman_res.pvalue
        spearman_ci_lo, spearman_ci_hi = compute_ci(spearman_r, len(est_cnt_lst))
        print(str(spearman_r) + "\t" + str(spearman_pval) + "\t" + str(spearman_ci_lo) +
              "\t" + str(spearman_ci_hi) + "\t")
    else:
        pearson_res = st.pearsonr(gt_cnt_lst, est_cnt_lst)
        pearson_r = pearson_res.statistic
        pearson_pval = pearson_res.pvalue
        pearson_ci = pearson_res.confidence_interval(confidence_level=0.95)
        print(str(pearson_r) + "\t" + str(pearson_pval) + "\t" + str(pearson_ci.low) +
              "\t" + str(pearson_ci.high) + "\t")


    # prepare density color
    xy = np.vstack([gt_cnt_lst, est_cnt_lst])
    z = st.gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = np.array(gt_cnt_lst)[idx], np.array(est_cnt_lst)[idx], z[idx]

    # draw scatter plot
    # plt.scatter(gt_cnt_lst, est_cnt_lst, alpha=0.5, s=10, color='black')
    plt.scatter(x, y, c=z, cmap='viridis', s=10)

    # Add color bar with logarithmic scaling
    cb = plt.colorbar()
    cb.set_label('Density')

    # Set color scale limits
    cb.set_ticks(np.linspace(z.min(), z.max(), 5))

    # linear regression
    slope, intercept = np.polyfit(gt_cnt_lst, est_cnt_lst, 1)
    plt.plot(gt_cnt_lst, np.polyval([slope, intercept], gt_cnt_lst), linestyle='-', color='red')

    if is_linear == 0:
        text_str = f'Spearman\'s $\\rho$ = {spearman_r:.3f}'
    else:
        text_str = f'Pearson\'s $r$ = {pearson_r:.3f}'

    plt.text(0.1, 0.9, text_str, transform=plt.gca().transAxes, backgroundcolor='none')

    # Add labels and title
    if is_linear == 0:
        plt.xlabel('log$_2$(expected read counts + 1)')
        plt.ylabel('log$_2$(estimated read counts + 1)')
    else:
        plt.xlabel('Expected read counts')
        plt.ylabel('Estimated read counts')

    # plt.xlim((-100, 5000))
    # plt.ylim((-100, 5000))
    plt.tight_layout()
    if is_svg == 1:
        plt.savefig(save_fn, format='svg')
    else:
        plt.savefig(save_fn)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
