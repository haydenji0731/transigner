#!/usr/bin/env python

import argparse
import sys
import pickle
import random
import os
from tqdm import tqdm
from subset import Tx, Gene, load_gtf


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-ip', '--index-prefix', type=str, help="", required=False, default=None)
    parser.add_argument('--gtf', type=str, help="", required=True)
    parser.add_argument('--gt-gtf', type=str, help="", required=False, default=None)
    parser.add_argument('--save', default=False, help="", required=False, action='store_true')
    parser.add_argument('-p', type=int, help="", required=True)
    parser.add_argument('-o', '--out-file', type=str, help="", required=True)
    parser.add_argument('--gt-only', default=False, help="", required=False, action='store_true')
    parser.add_argument('-prev', type=str, help="", required=False, default=None)
    parser.add_argument('--cont', default=False, help="", required=False, action='store_true')
    # TODO: implement w-genes mode for FLAIR2 benchmarking later...
    args = parser.parse_args()

    # perform argument checks
    if not 0 <= args.p <= 100:
        print("ARGUMENT ERROR: please provide a percentage argument that is > 1 and <= 100")
        sys.exit()

    if not args.gt_only and args.p == 0:
        print("percentage = 0 is only allowed in gt-only mode")

    if args.save and args.index_prefix is None:
        print("ARGUMENT ERROR: please provide an index path to save the transcriptome")
        sys.exit()

    if args.gt_only and args.gt_gtf is None:
        print("ARGUMENT ERROR: please provide the ground truth gtf file to use gt-only mode")
        sys.exit()

    if args.cont and args.prev is None:
        print("ARGUMENT ERROR: please provide the previously deleted transcripts list")
        sys.exit()

    skip = False

    if args.index_prefix is not None:
        fn = args.index_prefix + ".tx.index"
        if os.path.exists(fn):
            print("loading pre-built transcriptome index...")
            with open(fn, "rb") as f:
                tx_tbl = pickle.load(f)
            skip = True
    else:
        print("no pre-built transcriptome index detected... creating it")

    if not skip:
        tx_tbl = load_gtf(args.gtf)
        if args.save:
            print("saving the target transcriptome as an index...")
            fn = args.index_prefix + ".tx.index"
            with open(fn, "wb") as f:
                pickle.dump(tx_tbl, f)

    if args.cont:
        print("loading previously deleted transcripts list...")
        with open(args.prev, "rb") as f:
            prev_del = pickle.load(f)

    if args.gt_only:

        skip = False
        if args.index_prefix is not None:
            fn = args.index_prefix + ".gt.index"
            if os.path.exists(fn):
                print("loading pre-built ground truth transcriptome index...")
                with open(fn, "rb") as f:
                    gt_tx_tbl = pickle.load(f)
                skip = True
            else:
                if args.cont:
                    print("pre-built ground truth transcriptome index is required for this option")
                    sys.exit()
        if not skip:
            gt_tx_tbl = load_gtf(args.gt_gtf)
            if args.save:
                print("saving the ground truth transcriptome as an index...")
                fn = args.index_prefix + ".gt.index"
                with open(fn, "wb") as f:
                    pickle.dump(gt_tx_tbl, f)

        if args.p == 0:
            print('removing all ground truth transcripts')
            del_tid_lst = (gt_tx_tbl.keys())
        else:
            print("picking random indices to remove from the transcriptome")
            gt_tx_num = len(gt_tx_tbl)
            print("total # of ground truth transcripts: " + str(gt_tx_num))
            removal_perc = 100 - args.p
            print("percentage to remove: " + str(removal_perc) + "%")
            removal_num = int(gt_tx_num * removal_perc * 0.01)
            if args.cont:
                if len(prev_del) > removal_num:
                    print("the number of previously removed transcripts exceeds the number to be removed...")
                    print("please re-evaluate your threshold parameter")
                    sys.exit()
                elif len(prev_del) == removal_num:
                    print("the number of previously removed transcripts is equal to the number to be removed...")
                    print("please re-evaluate your threshold parameter")
                    sys.exit()
                removal_num -= len(prev_del)
                print("# of additional transcripts to remove on top: " + str(removal_num))
                rand_lst = random.sample(range(0, gt_tx_num), removal_num) + prev_del
            else:
                print("# of transcripts to remove: " + str(removal_num))
                rand_lst = random.sample(range(0, gt_tx_num), removal_num)

            if args.save:
                print("saving the removed transcripts list...")
                fn = args.index_prefix + ".removed.lst"
                with open(fn, "wb") as f:
                    pickle.dump(rand_lst, f)

            gt_tid_lst = list(gt_tx_tbl.keys())
            del_tid_lst = list()
            for i in tqdm(range(gt_tx_num)):
                tid = gt_tid_lst[i]
                if i in rand_lst:
                    del_tid_lst.append(tid)
    else:
        print("picking random indices to remove from the transcriptome")
        tx_num = len(tx_tbl)
        print("total # of transcripts: " + str(tx_num))
        removal_perc = 100 - args.p
        print("percentage to remove: " + str(removal_perc) + "%")
        removal_num = int(tx_num * removal_perc * 0.01)
        print("# of transcripts to remove: " + str(removal_num))
        if args.cont:
            if len(prev_del) > removal_num:
                print("the number of previously removed transcripts exceeds the number to be removed...")
                print("please re-evaluate your threshold parameter")
                sys.exit()
            elif len(prev_del) == removal_num:
                print("the number of previously removed transcripts is equal to the number to be removed...")
                print("please re-evaluate your threshold parameter")
                sys.exit()
            removal_num -= len(prev_del)
            rand_lst = random.sample(range(0, tx_num), removal_num) + prev_del
        else:
            rand_lst = random.sample(range(0, tx_num), removal_num)

        if args.save:
            print("saving the removed transcripts list...")
            fn = args.index_prefix + ".removed.lst"
            with open(fn, "wb") as f:
                pickle.dump(rand_lst, f)

    print("removing transcripts and writing the kept ones to a gtf file...")
    out_fh = open(args.out_file, 'w')

    if args.gt_only:
        tx_lst = list(tx_tbl.keys())
        for i in tqdm(range(len(tx_lst))):
            tid = tx_lst[i]
            if tid not in del_tid_lst:
                tx = tx_tbl[tid]
                out_fh.write(tx.raw + "\n")
                children = tx.children
                for child in children:
                    out_fh.write(child + "\n")
    else:
        tx_lst = list(tx_tbl.keys())
        for i in tqdm(range(len(tx_lst))):
            if i not in rand_lst:
                tid = tx_lst[i]
                tx = tx_tbl[tid]
                out_fh.write(tx.raw + "\n")
                children = tx.children
                for child in children:
                    out_fh.write(child + "\n")
    out_fh.close()


if __name__ == "__main__":
    main()










