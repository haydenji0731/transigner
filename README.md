# transigner

transigner is a program that assigns long RNA-seq reads to transcripts. Users can provide reads in FASTQ format and a set of transcript sequences in multi-FASTA format. [gffread](https://github.com/gpertea/gffread) can extract spliced out transcript sequences from annotation files in gtf/gff format.

## Installation

transigner is optimized for linux systems. We recommend that the users install the program in a conda environment.

```
// setting up conda env
$ conda create -n transigner python && conda activate transigner
$ conda install bioconda::minimap2 samtools

// easiest but without multithreading support
$ pip install transigner

// much faster (compiling from source)
$ pip install pysam numpy pandas
$ make CFLAGS="-fopenmp" && pip install .
```
Either clang or gcc must be available for compiling from source. transigner uses openmp for multithreading, so users should compile with the `-fopenmp` flag for optimal performance. By default, transigner installed through pip has multithreading disabled.

## Usage

transigner consists of three modules: align, prefilter, and em. easy.sh runs these modules at once:
```
easy.sh reads.fastq transcripts.fa out_dir // will execute align, pre, and em modules
```
There are static variables in this script that requires your attention:
```
aln_p # number of threads to use for minimap2 alignment
em_p # number of threas to use for EM iterations; only possible when openmp is available
data_type # input data type [ont_drna, ont_cdna, pacbio]
mode # transigner modes [default, psw, spiked]
```
We are still actively investigating the optimal parameter combinations for different data types, and the spiked mode is also under active development.

We provide a small set of query reads and target transcripts under `test_data/` directory. You'll get exactly 1 EM iteration. 

See below for the complete list of arguments:
```
Usage: transigner [MODULE] args [options]

*align module args / options:
      -q, --query
          query reads (FASTQ)
      -t, --target
          target transcripts (FASTA)
      -d, --out-dir
      -o, --out-file
      -n
          number of secondary alignments to retain
      -p, --threads
          number of threads to use for minimap2 alignment
*prefilter (pre) module args / options:
      -i, --in-file
          input alignment file (name sorted BAM file)
      -d, --out-dir
      --use-psw
          use position-specific weights when computing compatibility scores between reads and transcripts
*em (em) module args / options:
      -s, --scores
          compatibility scores (CSV)
      -d, --out-dir
      -u, --unmapped
          file containing unmapped read IDs
      -m, --tmap-file
          mappings between transcript IDs and integer indices (CSV)
      -c, --cvrg-thres
          EM convergence threshold defined in terms of cumulative read count change
      --push
          obtain hard 1-to-1 read assignments to transcripts
      --no-drop
          deactivate the drop feature that removes low scoring compatibility relationships
      -df, --drop-fac
          drop factor used to calculate the threshold for the drop
```

## Output

TODO

## Cite

> Ji, H. J. and M. Pertea (2024). "Enhancing transcriptome expression quantification through accurate assignment of long RNA sequencing reads with TranSigner." bioRxiv: 2024.2004.2013.589356. [doi:https://doi.org/10.1101/2024.04.13.589356]
