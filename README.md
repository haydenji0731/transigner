# TranSigner: assigning long RNA-seq reads to transcripts

TranSigner is a Python program that assigns long reads (compatible with both ONT and PacBio) to transcripts. This tool takes in a set of reads, and the assembled transcriptome. The transcriptome must be available in both fasta and gtf/gff formats. If you only have it in a gtf/gff format, you'll need to run [gffread](https://github.com/gpertea/gffread) as follows:
```
gffread -w transcripts.fa -g genome.fa transcripts.gtf
```
### Installation ###

You can install TranSigner with pip by running:
```
pip install transigner
```
or from source:
```
git clone https://github.com/haydenji0731/transigner transigner
cd transigner
python setup.py install
```
You also need [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/) installed. You can also use the precompiled binaries as long as they are available in your PATH. If you encounter trouble during installation, please open up a Github issue.

### Usage ###

TranSigner consists of three modules: align, prefilter, and em (short for expectation-maximization). You can check the arguments for each module by running:
```
transigner [module] -h
```
The `align` module takes a set of reads contained in a fastq file and align them to a transcriptome provided as a fasta file. See below:
```
transigner align -q reads.fastq -t transcripts.fa -d output_dir -o alignment.bam -p threads
```
Next, TranSigner processes the alignment results to compute compatibility scores between reads and transcripts and filter out alignments with questioning 5' and/or 3' end positions. The recommendedfilter thresholds differ by the read type:
```
transigner prefilter -a alignment.bam -t transcripts.fa -o output_dir --filter -tp -1 # noisy ONT direct RNA reads
transigner prefilter -a alignment.bam -t transcripts.fa -o output_dir --filter -tp -500 -fp -600 # ONT cDNA or PacBio IsoSeq
```
Finally, an expectation-maximization algorithm is run as follows:
```
transigner em -s output_dir/scores.tsv -i output_dir/ti.pkl -o output_dir --drop --use-score
```
By running all three modules, you'll obtain `abundances.tsv` and `assignments.tsv` files. See below for the header information for these files:
```
# abundances.tsv
transcript_id  read_count  relative_abundance
NR_024540.1	1.7149614376942495e-28	1.6757869165652874e-21

# assignments.tsv
read_id  (transcript_id_1, read_fraction_1)  (transcript_id_2, read_fraction_2)  (transcript_id_3, read_fraction_3) ...
57ebcba2-cde0-4096-b9b9-9bdc4306cb6c    (ENST00000559163, 6.640795584395568e-07)        (ENST00000559884, 0.9507564850356304)   (ENST00000354296, 0.04924285088481117)
```
You can also add the `--push` flag to obtain hard, 1-to-1 assignments between reads and transcripts.
### References ###
> Ji, H. J. and M. Pertea (2024). "Enhancing transcriptome expression quantification through accurate assignment of long RNA sequencing reads with TranSigner." bioRxiv: 2024.2004.2013.589356. [doi:https://doi.org/10.1101/2024.04.13.589356]

> Pertea, G., & Pertea, M. (2020). GFF utilities: GffRead and GffCompare [version 2; peer review: 3 approved].
> *F1000Research*, **9**:304. [doi:10.12688/f1000research.23297.2]

> Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *Gigascience*, > **10(2)**, giab008. [doi:10.1093/gigascience/giab008]

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
> *Bioinformatics*, **34**:3094-3100. [doi:10.1093/bioinformatics/bty191]

> Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
> *Bioinformatics*, **37**:4572-4574. [doi:10.1093/bioinformatics/btab705]

> Lianming Du, Qin Liu, Zhenxin Fan, Jie Tang, Xiuyue Zhang, Megan Price, Bisong Yue, Kelei Zhao. Pyfastx: a robust Python package for fast random access to sequences from plain and gzipped FASTA/Q files. Briefings in Bioinformatics, 2021, 22(4):bbaa368. 
