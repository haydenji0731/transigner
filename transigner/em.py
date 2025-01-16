from transigner.utils import *

def get_qt_sizes(fn):
    with open(fn, 'r') as fh:
        first_ln = fh.readline()
        temp = first_ln.strip()[1:].split(';')
        qsize = int(temp[0])
        tsize = int(temp[1])
    return qsize, tsize

def set_cvrg_thres(fn, dtype, tol, qsize, drop):
    ctr = 0
    with open(fn, 'r') as fh:
        for _ in fh: ctr += 1
    uperc = ctr / (ctr + qsize) * 100
    if dtype == 'ont_drna':
        t = 10 if uperc / tol[0] < 2 else 100
        if drop: t *= 10
    elif dtype == 'ont_cdna':
        # t = 10 if uperc / tol[1] < 2 else 100 # TODO: experiment with this
        t = 10
        if drop: t *= 10
    elif dtype == 'spiked':
        t = 10
    elif dtype == 'pacbio':
        t = 10
    else:
        print(tmessage("Unrecognized data type. Aborting...", Mtype.ERR))
        sys.exit(-1)
    return t
    
def main(args) -> None:
    check_dir(args.out_dir)
    param_fn = os.path.join(args.out_dir, "em_params.json")
    store_params(args, param_fn)

    qsize, tsize = get_qt_sizes(args.scores)

    if args.dtype == 'spiked' and not args.drop:
        print(tmessage("Running in spiked mode without --drop; We recommend using the drop feature", Mtype.PROG))

    if args.cvrg_thres == 'auto':
        cvrg_thres = set_cvrg_thres(args.unmapped, args.dtype, args.unmapped_tol, qsize, args.drop)
    else:
        cvrg_thres = args.cvrg_thres

    is_naive = 1 if args.naive else 0
    do_push = 1 if args.push else 0
    is_dev = 1 if args.dev else 0
    drop_fac = args.drop_fac if args.drop else 0
    relax_thres = args.relax_thres if args.relax else 0

    if is_naive:
        print(tmessage("NAIVE mode activated - scores calculated in the PRE module won't be used", Mtype.PROG))
    if do_push:
        print(tmessage("PUSH flag specified - hard 1-to-1 assignments will be outputted", Mtype.PROG))

    cmd = (
        f'em -s {args.scores} '
        f'-m {args.tmap_file} '
        f'-qs {qsize} '
        f'-ts {tsize} '
        f'-p {args.threads} '
        f'--naive {is_naive} '
        f'-n {args.num_iter} '
        f'-c {cvrg_thres} '
        f'-o {args.out_dir} '
        f'--push {do_push} '
        f'--dev {is_dev} '
        f'-r {relax_thres} '
        f'-df {drop_fac}'
    )
    print(cmd)
    call(cmd, shell=True)
