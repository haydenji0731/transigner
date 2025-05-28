## Getting Started

```
# compile from source
git clone --recurse-submodules https://github.com/haydenji0731/transigner
cd transigner && make release

./transigner aln.bam # alignment-based mode
./transigner -t txp.fa reads.fq # read-based mode; uses minimap2

# help msg contains command line options
./transigner -h
```

## Table of Contents

- Installation
- Usage
  - Alignment-based mode
  - Read-based mode

## Installation

transigner is optimized for x86-64 CPU and uses OpenMP for parallelization.


## Usage
### Alignment-based mode
