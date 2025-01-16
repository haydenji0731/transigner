#!/usr/bin/env bash

reads=$1
transcripts=$2
out_dir=$3

aln_p=24 # number of threads to use for minimap2 alignment
em_p=4 # number of threads to use for EM iterations
data_type="ont_drna" # input data type
mode="default"

set -x

transigner align -q $reads -t $transcripts -d $out_dir \
   -o "aligned.bam" -p $aln_p
if [ $mode == "default" ]; then
    transigner pre -i "${out_dir}/aligned.bam" -d $out_dir 
    transigner em -s "${out_dir}/scores.csv" -d $out_dir \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype $data_type -p $em_p
elif [ $mode == "psw"]; then
    transigner pre -i "${out_dir}/aligned.bam" -d $out_dir --use-psw
    transigner em -s "${out_dir}/scores.csv" -d $out_dir \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype $data_type -p $em_p --no-drop
elif [ $mode == "spiked" ]; then
    transigner pre -i "${out_dir}/aligned.bam" -d $out_dir --spiked
    transigner em -s "${out_dir}/scores.csv" -d $out_dir \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype $data_type -p $em_p
else
    echo "unrecognized mode"
fi