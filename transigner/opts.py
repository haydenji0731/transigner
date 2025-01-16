index_opt_d = {
    '-k': "15",
    '-w': "10",
    '-H': False,
    '-I': '4G',
    '--idx-no-seq': False,
    '-d': None,
    '--alt': None,
    '--alt-drop': "0.15"
}

map_opt_d = {
    '-f': '0.0002',
    '-U': '10, 1000000',
    '--q-occ-frac': '0.01',
    '-e': '500',
    '-g': '10k',
    '-r': '500, 20k',
    '-n': '3',
    '-m': '40',
    '-D': False,
    '-P': False,
    '--dual': 'yes',
    '-X': False,
    '-p': '0.8',
    '-N': '5',
    '-G': '200k',
    '-F': '800',
    '-M': '0.5',
    '--rmq': 'no',
    '--hard-mask-level': False,
    'mask-len': 'inf',
    '--max-chain-skip': '25',
    '--max-chain-iter': '5000',
    '--chain-gap-scale': '1.0',
    '--no-long-join': False,
    '--splice': False,
    '--sr': False,
    '--split-prefix': None,
    '--frag': 'no',
    '--for-only': False,
    '--rev-only': False,
    '--heap-sort': 'no',
    '--no-pairing': False
}

baln_opt_d = {
    '-A': '2',
    '-B': '4',
    '-O': '4, 24',
    '-E': '2, 1',
    '-C': '0',
    '-z': '400, 200',
    '-s': '40',
    '-u': 'n',
    '--end-bonus': '0',
    '--score-N': '1',
    '--splice-flank': 'no',
    '--junc-bed': None,
    '--junc-bonus': '9',
    '--end-seed-pen': '6',
    '--no-end-flt': False,
    '--cap-sw-mem': '100m',
    '--cap-kalloc': '0'
}

io_opt_d = {
    '-a': False,
    '-o': None,
    '-Q': False,
    '-L': False,
    '-R': None,
    '-y': False,
    '-c': False,
    '--cs': False,
    '--MD': False,
    '--eqx': False,
    '-Y': False,
    '--seed': '11',
    '-t': '3',
    '-2': False,
    '-K': '500M',
    '--secondary':'yes',
    '--max-qlen': None,
    '--paf-no-hit': False,
    '--sam-hit-only': False
}

presets = [
    'map-ont',
    'map-hifi',
    'map-pb',
    'asm5',
    'asm10',
    'asm20',
    'splice',
    'splice:hq',
    'sr',
    'ava-pb',
    'ava-ont'
]

def opt_dict2str(opt_d):
    opt_s = []
    for k in opt_d:
        if opt_d[k] is True:
            opt_s.append(f'{k}')
            continue
        elif opt_d[k] is False or opt_d[k] is None:
            continue
        else:
            opt_s.append(f'{k} {opt_d[k]}')
    opt_s = ' '.join(opt_s)
    return f'{opt_s} '