## Getting Started

```
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
- Output

## Installation

transigner is optimized for x86-64 CPU and uses OpenMP for parallelization.


## Usage
### Alignment-based mode

If you already have an alignment BAM file, simply provide it as follows:

```
./transigner aln.bam
```

### Read-based mode

If you have a FASTQ file of reads to process, prepare a FASTA file containing transcript sequences to quantify / assign to. [gffread](https://github.com/gpertea/gffread) can be useful here:

```
# compile gffread from source (preferred) or install from bioconda
gffread -w txp.fa -g genome.fa txp.gff/gtf
```

transigner calls minimap2 to align reads to transcripts so install minimap2 if you haven't already. Then call:

```
./transigner -t txp.fa reads.fastq
```

## Output
