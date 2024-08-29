#!/usr/bin/env python

import os
import json
from datetime import datetime
import pickle
import sys
from transigner.commons import RED, GREEN, RESET

def load_scores(fn, ti_map):
    qi_map = dict()
    ti_set = set()
    first = True
    qi = -1
    with open(fn, 'r') as f:
        for ln in f:
            if first:
                qi_size = int(ln.replace('#', ''))
                cmpt_mat = [dict() for _ in range(qi_size)]
                first = False
                print(datetime.now(), f"{GREEN}PROGRESS{RESET} number of reads:", qi_size)
                continue

            fields = ln.strip(). split("\t")
            for i in range(len(fields)):
                if i == 0:
                    qi += 1
                    qname = fields[i]
                    qi_map[qi] = qname # different from prefilter
                else:
                    tname = fields[i].split(',')[0].replace('(', '')
                    score = float(fields[i].split(",")[1].strip().replace(')', ''))
                    ti = ti_map[tname]
                    cmpt_mat[qi][ti] = score
                    ti_set.add(ti)
    # sanity check
    assert qi_size == len(qi_map)
    return cmpt_mat, qi_map, qi_size, ti_set


def init(qi_size, ti_set, ti_size, pre_init, p_est=None):
    rho = dict()
    alpha = [dict() for _ in range(qi_size)]
    if pre_init:
        sub_sum = 0
        for ti in ti_set:
            sub_sum += p_est[ti]
        for ti in ti_set:
            rho[ti] = p_est[ti] / sub_sum
    else:
        for ti in ti_set:
            rho[ti] = 1 / ti_size
    return rho, alpha

def step_e(alpha, rho, qi_size, cmpt_mat, use_score):
    for qi in range(qi_size):
        alpha_qi = dict()
        loc_ti_set = cmpt_mat[qi]
        for ti in loc_ti_set:
            score = cmpt_mat[qi][ti]
            rho_ti = rho[ti]
            if use_score:
                alpha_qi[ti] = score * rho_ti
            else:
                alpha_qi[ti] = rho_ti
        z = sum(alpha_qi.values())
        # TODO: is this correct?
        if z == 0:
            for ti in loc_ti_set:
                alpha[qi][ti] = 0
        else:
            for ti in loc_ti_set:
                alpha[qi][ti] = alpha_qi[ti] / z

def step_m(alpha, ti_set, qi_size):
    new_rho = dict()
    for ti in ti_set:
        new_rho[ti] = 0
    for qi in range(qi_size):
        for ti in alpha[qi]:
            new_rho[ti] += alpha[qi][ti]
    for ti in ti_set:
        new_rho[ti] /= qi_size # a read == a transcript
    return new_rho

def has_converged(old_rho, new_rho, thres):
    delt_rho = 0
    converged = False
    for ti in old_rho:
        # delt_rho += abs(old_rho[ti] - new_rho[ti]
        delt_rho = max(abs(old_rho[ti] - new_rho[ti]), delt_rho)
    if delt_rho < thres:
        converged = True
    return delt_rho, converged
        
def drop_scores(cmpt_mat, alpha, qi_size, df):
    for qi in range(qi_size):
        n_qi = sum(1 for s in cmpt_mat[qi].values() if s > 0)
        if n_qi <= 1:
            continue
        sigma_qi = 1 / n_qi + (1 / n_qi * df)
        max_alpha = max(cmpt_mat[qi].values())
        max_tis = [k for k, v in cmpt_mat[qi].items() if v == max_alpha]
        ctr = 0
        for ti in cmpt_mat[qi]:
            rf = alpha[qi][ti]
            if rf < sigma_qi:
                cmpt_mat[qi][ti] = 1e-2 # no fraction of qi assigned to ti
            else:
                ctr += 1
        if ctr == 0:
            for ti in max_tis:
                cmpt_mat[qi][ti] = max_alpha

def relax(rho, ti_set, qi_size):
    new_rho = dict()
    for ti in ti_set:
        if rho[ti] < (1 / qi_size):
            # print(rho[ti])
            new_rho[ti] = 0
        else:
            new_rho[ti] = rho[ti]
    return new_rho

def load_pre_est(fn, ti_map):
    p_est = [0] * len(ti_map)
    with open(fn, 'r') as f:
        for ln in f:
            fields = ln.strip().split("\t")
            tname = fields[0]
            try:
                assert tname in ti_map
            except:
                print(datetime.now(), f"{RED}ERROR{RESET} unrecognized transcript \
                      specified in the pre-estimate")
                sys.exit(-1)
            ti = ti_map[tname]
            est_ti = float(fields[1])
            p_est[ti] = est_ti
    return p_est


def main(args):
    if args.pre_init and args.estimate is None:
        print(datetime.now(), f"{RED}ERROR{RESET} please provide a pre-estimate file")
        sys.exit(-1)

    if args.use_score and args.scores is None:
        print(datetime.now(), f"{RED}ERROR{RESET} please provide a scores file")
        sys.exit(-1)

    if args.naive and (args.use_score or args.drop):
        print(datetime.now(), f"{RED}ERROR{RESET} can't use score and/or drop in naive mode")
        sys.exit(-1)
    
    cmd_fn = os.path.join(args.out_dir, "em_cmd_info.json")
    with open(cmd_fn, 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading target index")
    with open(args.index, 'rb') as f:
        ti_map = pickle.load(f)
    
    if args.pre_init:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading pre-estimates")
        p_est = load_pre_est(args.estimate, ti_map)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} loading compatibility scores")
    cmpt_mat, qi_map, qi_size, ti_set = load_scores(args.scores, ti_map)
    ti_size = len(ti_set)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} initializing")
    if args.pre_init:
        rho, alpha = init(qi_size, ti_set, ti_size, True, p_est)
    else:
        rho, alpha = init(qi_size, ti_set, ti_size, False)

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} starting EM")
    num_iter = 1
    rho_converged = False
    while num_iter <= args.max_iter:
        step_e(alpha, rho, qi_size, cmpt_mat, args.use_score)
        if args.drop:
            if num_iter == 1:
                drop_scores(cmpt_mat, alpha, qi_size, args.drop_fac)
                step_e(alpha, rho, qi_size, cmpt_mat, args.use_score)
        new_rho = step_m(alpha, ti_set, qi_size)
        delta_rho, converged = has_converged(rho, new_rho, args.rho_thres)
        if args.verbose:
            print("iteration", num_iter, "cumulative rho delta:", delta_rho)
        rho = new_rho
        if args.relax:
            new_rho = relax(rho, ti_set, qi_size)
            rho = new_rho
        if converged:
            print(datetime.now(), f"{GREEN}PROGRESS{RESET} converged")
            rho_converged = True
            break
        num_iter += 1
    if not rho_converged:
        print(datetime.now(), f"{GREEN}PROGRESS{RESET} max iter reached")
    
    summary_fn = os.path.join(args.out_dir, "em_summary.txt")
    with open(summary_fn, 'w') as f:
        f.write(f"number of EM iterations: {num_iter - 1}\n")
        f.write(f"last cumulative rho delta: {delta_rho}\n")
    
    alpha_fn = os.path.join(args.out_dir, "assignments.tsv")
    un_asgn_fn = os.path.join(args.out_dir, "unassigned.txt")
    un_asgn_lines = []
    alpha_lines = []

    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out assignments")
    tnames = sorted(ti_map, key=lambda k: ti_map[k])
    for qi in tqdm(range(qi_size)):
        qname = qi_map[qi]
        if sum(alpha[qi].values()) == 0:
            print(datetime.now(), f"{RED}WARNING{RESET} read {qname} unassigned")
            un_asgn_lines.append(f"{qname}\n")
            continue
        line = [qname]
        for ti, alpha_ti in alpha[qi].items():
            if alpha_ti > 0:
                tname = tnames[ti]
                line.append(f"\t({tname}, {alpha_ti})")
        alpha_lines.append("\t".join(line) + "\n")
    
    with open(alpha_fn, 'w') as alpha_fh:
        alpha_fh.writelines(alpha_lines)
    
    with open(un_asgn_fn, 'w') as un_asgn_fh:
        un_asgn_fh.writelines(un_asgn_lines)

    if args.push:
        pushed_alpha_fn = os.path.join(args.out_dir, "hard_assignments.tsv")
        pushed_alpha_fh = open(pushed_alpha_fn, 'w')
        for qi in range(qi_size):
            qname = qi_map[qi]
            if sum(alpha[qi].values()) == 0:
                continue
            pushed_alpha_fh.write(qname)
            alpha_qi = alpha[qi]
            max_alpha = max(alpha_qi.values())
            max_tis = {key for key, value in alpha_qi.items() if value == max_alpha}
            if len(max_tis) < 1:
                print(datetime.now(), f"{RED}WARNING{RESET} a tie detected for read {qname}")
            for ti in max_tis:
                tname = tnames[ti]
                pushed_alpha_fh.write("\t" + tname)
            pushed_alpha_fh.write("\n")
        pushed_alpha_fh.close()
    
    print(datetime.now(), f"{GREEN}PROGRESS{RESET} writing out abundances")
    rho_fn = os.path.join(args.out_dir, "abundances.tsv")
    rho_lines = []
    tnames_set = set(tnames)
    for ti in ti_set:
        tname = tnames[ti]
        rho_ti = rho[ti]
        rc_ti = rho_ti * qi_size
        rho_lines.append(f"{tname}\t{rho_ti}\t{rc_ti}\n")

    for tname in ti_map:
        if tname not in tnames_set:
            rho_lines.append(f"{tname}\t0.0\t0\n")
    
    with open(rho_fn, 'w') as rho_fh:
        rho_fh.writelines(rho_lines)
        
if __name__ == "__main__":
    main()