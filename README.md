# TranSigner: assigning long RNA-seq reads to assembled transcripts

TranSigner is a Python program that assigns long ONT reads to the assembled transcripts. This tool takes in a set of reads, and the assembled transcriptome. The transcriptome must be available in both Fasta and GTF/GFF formats. If you only have it in a GTF/GFF format, you'll need to run [GffRead](https://github.com/gpertea/gffread) as follows:

```
$ gffread -w assembled_txs.fa -g genome.fa assembled_txs.gtf
```
Make sure to use the same version of genome that you used to assembled transcripts. 

### Installation ###

TranSigner is divided up into 3 stages: align, prefilter, and em (short for expectation-maximization). These steps internally calls [minimap2](https://github.com/lh3/minimap2), [samtools](http://www.htslib.org/), and [jf_aligner](https://github.com/alekseyzimin/jf_aligner). It is required that you install these tools and have it available either in a conda environment and/or in your local path at the time of running TranSigner. I'm suggesting one possible way to install these softwares in your system, but unfortunately, there might be some errors specific to the Apple's ARM silicon chip, so expect to do some debugging if you are a Mac user. In the future, I plan to provide an ARM processor Q&A section in this README.

Assuming you have [miniconda](https://docs.conda.io/projects/miniconda/en/latest/) installed in your system, I'm including one possible way to install these tools in your conda environment:

1. Installing minimap2
```
$ curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
$ ./minimap2-2.26_x64-linux/minimap2
$ cd minimap2-2.26_x64-linux
$ export PATH=$PATH:`pwd`
$ minimap2 -V
```
If this doesn't return the minimap2 manual page, then it might be that (1) the precompiled binary is not compatible with your system or (2) your export command failed for some unknown reason.

2. Installing samtools
```
$ conda install -c bioconda samtools=1.9
$ samtools
```
If this doens't return the samtools manual page, then your samtools installation through conda package manager must have failed. You can attempt building it form the source following instructions [here](http://www.htslib.org/download/). I would recommend that you use version 1.9 instead of 1.18.
  
4. Installing jf_aligner (OPTIONAL)
jf_aligner requires a [Boost](https://www.boost.org/) installation. If you already have Boost installed in your system, you can run
```
$ export BOOST_ROOT=/path/to/boost/
```
before running the jf_aligner installation script. Otherwise, you can also request that Boost be installed on the fly as part of this script by running:
```
$ export BOOST_ROOT=install
```
Now, you can download the jf_aligner package, unzip it, and run the installation script as follows:
```
$ wget https://github.com/alekseyzimin/jf_aligner/releases/download/v1.0.1/jf_aligner-1.0.1.tar.gz
$ tar -xf jf_aligner-1.0.1.tar.gz
$ cd jf_aligner-1.0.0
$ ./install.sh
$ jf_aligner -h
```

Lastly, we need to install few Python packages required to run TranSigner. For this, there's a `requirements.txt` file included in this git repo that you can use to install these packages all at once. Note that if you've installed all above tools in your local environment without conda, you can also create a Python [virtual environment](https://docs.python.org/3/library/venv.html) and run the same exact commands as following:

```
$ git clone https://github.com/haydenji0731/TranSigner.git
$ cd TranSigner
$ pip install -r requirements.txt
```

### How to Run ###

A single TranSigner run consists of running 3 separate Python scripts: `align.py`, `prefilter.py`, and `em.py` in this order. The expected runtime for completing all these modules ~1h depending on the read set size and number of threads used. At the end, the program will output 2 files: `abundance.tsv` and `assignment.tsv`. 

In theory, we can skip the jf_aligner call in prefilter module to expedite our TranSigner runs, so if you are attempting to process multiple files and time is a big concern, then please open a issue asking that I make the jf_aligner call be optional. 

Now, we can run the `align` module as follows:

```
$ python ./align.py -i reads.fastq -ref-fa assembled_txs.fasta -o aln_output_dir -op output_prefix -t num_threads
```
If you are expecting a high number of isoforms at some gene loci, then you can increase the `-sN` optional parameter. Currently, the default value for this parameter is set of 181, which was calibrated based on the most recent RefSeq hg38 human genome annotation.

Next, run the following command to start the `prefilter` module:

```
# default - without jf_aligner call
$ python ./prefilter.py -query reads.fastq -target assembled_txs.fasta -aln ${aln_output_dir}/aln.bam -gtf assembled_txs.gtf -t num_threads -tmp ${pref_output_dir}/tmp -o pref_output_dir --skip-jf --filter
$ python ./prefilter.py -query reads.fastq -target assembled_txs.fasta -aln ${aln_output_dir}/aln.bam -gtf assembled_txs.gtf -t num_threads -tmp ${pref_output_dir}/tmp -o pref_output_dir --score-ratio X Y --filter 
```
`--score-ratio` indicates how much weight is assigned to minimap vs. jf_aligner. I recommend that 99% of the weight goes to the minimap output and 1% goes to jf_aligner (i.e., --score-ratio 0.99 0.01). This achieves a small increase in our read assignment performance. 

Also, you can use `--five-prime` and `--three-prime` arguments to specify how lenient you want to be with the distance between the 5' end of a transcript and mapping start position of a read, when filtering out non-primary (i.e., secondary, supplementary alignment). By default, I allow upto 650bp difference with the primary alignment's end distances on the 5' end. 

We'll proceed to run the last module. Before running it, we need to know the size of your compatibility matrix constructed during the prefilter module by running:
```
$ wc -l ${pref_output_dir}/cmpt_tbl.tsv
8920552 (this number will be different for you)
```
Take note of this number since you'll need it for running the actual `em` module. Afterwards, you can run `em` by running:
```
$ ./em.py -cmpt ${pref_output_dir}/cmpt_tbl.tsv -scores ${pref_output_dir}/scores.tsv -i ${pref_output_dir}/transcriptome.index --use-score --drop -o $em_output_dir -l 8920552 (this number will be different for you)
```
### References ###

> Zimin, A. V., Marçais, G., Puiu, D., Roberts, M., Salzberg, S. L., & Yorke, J. A. (2013). The MaSuRCA genome assembler.
> *Bioinformatics*, **29(21)**:2669-2677. [doi:10.1093/bioinformatics/btt476]

> Pertea, G., & Pertea, M. (2020). GFF utilities: GffRead and GffCompare [version 2; peer review: 3 approved].
> *F1000Research*, **9**:304. [doi:10.12688/f1000research.23297.2]

> Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *Gigascience*, > **10(2)**, giab008. [doi:10.1093/gigascience/giab008]

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
> *Bioinformatics*, **34**:3094-3100. [doi:10.1093/bioinformatics/bty191]

> Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
> *Bioinformatics*, **37**:4572-4574. [doi:10.1093/bioinformatics/btab705]
