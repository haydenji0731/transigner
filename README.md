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
- Arguments

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

transigner calls minimap2 to align reads to transcripts so install it if you haven't already. Then call:

```
./transigner -t txp.fa reads.fastq [-dtype ont or pacbio]
```

Specifying `-dtype` changes the minimap2 preset used during alignment.

## Output

- `alignments.bam`: only generated in read-based mode. minimap2-generated bam file.

- `assignments.tsv`: read-to-transcript assignment file with the following format:

```
read_id<TAB>transcript_id,fraction_assigned<TAB>transcript_id,fraction_assigned...
```
Each line represents a single read and lists the transcripts it was assigned to along with the fraction of the read assigned to each transcript. Fractions sum to 1 (or close to 1) across all assigned transcripts for that read.

- `abundances.csv`: transcript-level read counts with the following format:

```
transcript_id,read_count
```

transcripts with no reads aligned to get assigned with a count of 0.

## Additional arguments

```
Usage:
  transigner [OPTION...] input_file

  -o, --out-dir arg    output directory (default: .)
  -t, --tx arg         fasta file containing transcript sequences
  -a, --max-aln arg    maximum number of alignments per read (default: 181)
  -n, --max-iter arg   number of max EM iterations (default: 1000)
  -e, --epsilon arg    EM convergence threshold ε (default: 10.0)
  -c, --constant arg   τ constant in exponential decay fn (default: 5.0)
      --sample         assign by sampling (default: false)
      --use-psw        use positional weights (default: false)
  -r, --min-rho arg    zero out txes with read counts below this threshold 
                       (default: 0.1)
      --keep-low       keep low-prob alignments (default: false)
  -k, --min-pk arg     constant to computed min alignment prob (default: 
                       0.1)
  -w, --base-w arg     base weight for an alignment (default: 0.0)
  -d, --data-type arg  long read RNA-seq data type (default: ont [ont or pacbio accepted])
  -p, --threads arg    number of threads (default: 1)
  -v, --verbose        enable verbose output for more detailed logging 
                       (default: false)
  -h, --help           usage
```
