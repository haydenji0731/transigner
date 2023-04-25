# TranSigner: assigning long reads to transcripts

A brief Python program to assigns long ONT reads to the assembled transcripts. This tool takes reads and reference transcriptome (FASTA and GTF) as required inputs. 

### Installation ###

TranSigner is divided up into 2 stages: align and em. The align module internally calls [minimap2](https://github.com/lh3/minimap2) and [stringtie2](https://github.com/gpertea/stringtie). It is required that you install these tools and have it available in your local PATH. I intend to update the requirements.txt to allow users to create a conda environment to download all dependencies at once. For now, the file only contains Python package dependencies and thus simply creating a Python virtual environment would suffice. To do so run the following commands:

```
$ git clone https://github.com/haydenji0731/TranSigner.git
$ cd TranSigner
$ python3 -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```

### How to Run ###

As mentioned before, TranSigner first aligns the input reads to the reference transcriptome and then assigns each read to a transcript through an iterative EM execution. The expected runtime for completing both modules is 1~2 hours. At the end, the program will output 2 files: ```abundance.tsv``` and ```assignment.tsv```. Run the following command to begin the alignment stage:

```
$ python ./align.py -i reads.fastq -ref-fa transcriptome.fasta -o output_dir
```

Assuming that you've already assembled a transcriptome with your input reads, the corresponding FASTA file can be obtained using [gffread](https://github.com/gpertea/gffread). Above command also shows the required positional arguments. Optional arguments include: number of threads (```-t```) and output prefix (```-op```). To learn more about these optional parameters, simply run:

```
$ python ./align.py -h
```

Next, run the following command to begin the EM stage:
```
$ python ./em.py -i reads.fastq -ref-gtf transcriptome.gtf -o output_dir
```
Above command also shows the required positional arguments. Optional arguments include: min TPM change threshold for convergence (```-thres```), maximum number of iterations (```-max-iter```), and output prefix (```-op```). To learn more about these optional parameters, simply run:

```
$ python ./em.py -h
```

Finally, there are some auxiliary Python scripts for computing mean absolute TPM loss or precision. See ```./scripts``` folder for more information. 
