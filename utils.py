#!/usr/bin/env python

from tqdm import tqdm
from operator import itemgetter
import sys


def build_tbls_paf(reads_lst, reads_tbl, tx_tbl, padding=1):
    reads_num = len(reads_lst)
    cover_tbl = [dict() for _ in range(reads_num)]
    cmpt_tbl = [set() for _ in range(reads_num)]

    for i in tqdm(range(reads_num)):
        qname = reads_lst[i]
        recs = reads_tbl[qname]
        cover_lst = list()
        for rec in recs:
            fields = rec.split(" ")
            tname = fields[14]
            t_idx = tx_tbl[tname]
            cmpt_tbl[i].add(t_idx)

            r_st = int(fields[0])
            r_en = int(fields[1])
            q_st = int(fields[2])
            q_en = int(fields[3])
            r_cover = r_en - r_st
            q_cover = q_en - q_st
            r2q_cover_ratio = r_cover / q_cover
            r2q_cover_delt = -1 * abs(r2q_cover_ratio - 1)
            cover_lst.append((t_idx, r2q_cover_delt))

        shift = abs(min(cover_lst, key=itemgetter(1))[1]) + padding
        cover_tmp_lst = [(i, x + shift) for (i, x) in cover_lst]
        cover_max = max(cover_tmp_lst, key=itemgetter(1))[1]
        cover_norm_lst = [(i, x / cover_max) for (i, x) in cover_tmp_lst]

        for ti, x in cover_norm_lst:
            cover_tbl[i][ti] = x

    return cover_tbl, cmpt_tbl


def merge_cmpt_tbl(cmpt_tbl_1, reads_lst_1, cmpt_tbl_2, reads_lst_2, ms_tbl, nmer_tbl, stretch_tbl, offset_tbl,
                   cover_tbl, padding=0.01):
    # TODO: consider pursuing a full UNION approach
    # TODO: test out different padding parameters
    print("preparing for the merge...")
    reads_num_1 = len(reads_lst_1)
    reads_num_2 = len(reads_lst_2)
    reads_tbl_2 = dict()
    for i in tqdm(range(reads_num_2)):
        r = reads_lst_2[i]
        reads_tbl_2[r] = i

    mg_reads_lst = list()
    print("merging read lists first...")
    for i in tqdm(range(reads_num_1)):
        r = reads_lst_1[i]
        if r in reads_tbl_2:
            j = reads_tbl_2[r]
            mg_reads_lst.append((r, i, j))
        else:
            mg_reads_lst.append((r, i, -1))

    print("merging compatibility matrices...")
    mg_reads_num = len(mg_reads_lst)
    mg_cmpt_tbl = [set() for _ in range(mg_reads_num)]
    mg_scores_mat = [dict() for _ in range(mg_reads_num)]
    # score = ms * nmer * stretch * offset * cover
    for k in tqdm(range(mg_reads_num)):
        qname, i, j = mg_reads_lst[k]
        if j == -1:
            ti_set = cmpt_tbl_1[i]
            mg_cmpt_tbl[k] = ti_set
            for ti in ti_set:
                ms = ms_tbl[i][ti]
                score = ms * padding
                # score = ms * (padding ** 4)
                mg_scores_mat[i][ti] = score
        else:
            ti_set_1 = cmpt_tbl_1[i]
            ti_set_2 = cmpt_tbl_2[j]
            mg_ti_set = ti_set_1.union(ti_set_2)
            mg_cmpt_tbl[k] = mg_ti_set
            for ti in mg_ti_set:
                is_set_1 = ti in ti_set_1
                is_set_2 = ti in ti_set_2
                if is_set_1 and is_set_2:
                    ms = ms_tbl[i][ti]
                    nmer = nmer_tbl[j][ti]
                    stretch = stretch_tbl[j][ti]
                    offset = offset_tbl[j][ti]
                    cover = cover_tbl[j][ti]
                    score = ms * nmer * stretch * offset * cover
                elif is_set_1:
                    ms = ms_tbl[i][ti]
                    score = ms * padding
                elif is_set_2:
                    nmer = nmer_tbl[j][ti]
                    stretch = stretch_tbl[j][ti]
                    offset = offset_tbl[j][ti]
                    cover = cover_tbl[j][ti]
                    score = padding * nmer * stretch * offset * cover
                else:
                    print("FATAL ERROR: impossible edge case reached.")
                    sys.exit()
                if score < 0:
                    print("FATAL ERROR: score cannot be negative.")
                    sys.exit()
                mg_scores_mat[i][ti] = score

    return mg_cmpt_tbl, mg_scores_mat, mg_reads_lst


def build_score_tbl(cmpt_tbl_1, reads_lst_1, cmpt_tbl_2, reads_lst_2, ms_tbl, offset_tbl, ratio):
    print("preparing to build the score matrix")
    reads_num_1 = len(reads_lst_1)
    reads_num_2 = len(reads_lst_2)
    reads_tbl_2 = dict()
    score_tbl = [dict() for _ in range(reads_num_1)]

    for i in tqdm(range(reads_num_2)):
        r = reads_lst_2[i]
        reads_tbl_2[r] = i

    print("building the score matrix...")
    for i in tqdm(range(reads_num_1)):
        qname = reads_lst_1[i]
        ti_set_1 = cmpt_tbl_1[i]
        if qname in reads_tbl_2:
            j = reads_tbl_2[qname]
            ti_set_2 = cmpt_tbl_2[j]
            for ti in ti_set_1:
                if ti in ti_set_2:
                    score = ms_tbl[i][ti] * ratio[0] + offset_tbl[j][ti] * ratio[1]
                else:
                    score = ms_tbl[i][ti]
                score_tbl[i][ti] = score
        else:
            for ti in ti_set_1:
                score = ms_tbl[i][ti]
                score_tbl[i][ti] = score

    return score_tbl


