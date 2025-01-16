from transigner.utils import *
from transigner.opts import *

def process_opts(in_s, opt_d):
    if in_s == "": # no option specified
        return ""
    fields = in_s.replace('"', '').strip().split()
    i = 0
    new_opt_d = dict()
    while i < len(fields):
        f = fields[i]
        if f in opt_d:
            if opt_d[f] is False:
                new_opt_d[f] = True
            else:
                if i + 1 >= len(fields):
                    print(tmessage("Couldn't process options", Mtype.ERR))
                    sys.exit(-1)
                new_opt_d[f] = fields[i+1]
        i += 1
    return(opt_dict2str(new_opt_d))

def main(args):
    param_fn = os.path.join(args.out_dir, "align_params.json")
    store_params(args, param_fn)
    ofn = os.path.join(args.out_dir, args.out_file)
    if args.dev:
        index_opt_s = process_opts(args.index_opts, index_opt_d)
        map_opt_s = process_opts(args.map_opts, map_opt_d)
        baln_opt_s = process_opts(args.base_aln_opts, baln_opt_d)
        cmd = f'minimap2 -ax {index_opt_s}{map_opt_s}{baln_opt_s}{args.query} {args.target} {args.threads}'
    else:
        if args.preset not in presets:
            print(tmessage("Unknown mm2 preset", Mtype.ERR))
            sys.exit(-1)
        cmd = f'minimap2 -ax {args.preset} -N {args.n} -t {args.threads} {args.target} {args.query}'
    cmd += f' | samtools view -b -o {ofn} -@ {args.threads}'
    print(tmessage("Alignment started", Mtype.PROG))
    print(cmd)
    call(cmd, shell=True)