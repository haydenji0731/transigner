#/bin/bash

in_fq=$1
op=$2
threads=$3
resultDir="/ccb/salz8-3/hji20/ela/transigner/results"
dataDir="/ccb/salz8-3/hji20/ela/transigner/data"
readDir="/ccb/salz8-3/hji20/ela/transigner/results/nanosim_reads"

echo "input fastq file: " $in_fq
minimap2 -ax splice -t $threads -k14 "${dataDir}/hg38.p13.fa" "${readDir}/${in_fq}" > "${resultDir}/minimap2/${op}.sam"
samtools view -Sb -o "${resultDir}/minimap2/${op}.bam" "${resultDir}/minimap2/${op}.sam" -@ $threads
samtools sort -o "${resultDir}/minimap2/${op}.sorted.bam" "${resultDir}/minimap2/${op}.bam" -@ $threads
samtools index "${resultDir}/minimap2/${op}.sorted.bam" -@ $threads
stringtie -L -G "${dataDir}/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation_no_genes.gtf" -o "${resultDir}/stringtie2/${op}.gtf" "${resultDir}/minimap2/${op}.sorted.bam" -p $threads -l ${op}