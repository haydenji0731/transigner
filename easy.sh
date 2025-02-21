#!/usr/bin/env bash

# Set flag defaults
a_flag=0
e_flag=0
d_flag=0
m_flag=0

##########################################################
#                        HELP                            #
##########################################################
function Help() 
{   
# Show help 
    echo "transigner start-to-end in one command - use correct order of positional arguments.
    Usage: 
        easy.sh -a aligner_threads -e em_threads -d data_type -m mode <reads> <transcripts> <outdir> 
        easy.sh -h 
        easy.sh -v  
    
    Arguments:
        reads           input reads
        transcripts     reference transcripts
        outdir          output directory
        -a              minimap2 alignment threads
        -e              EM iterations threads (requires openmp during compilation)
        -d              data type [ont_drna, ont_cdna, pacbio]
        -m              transigner mode [default, psw, spiked]
    
    Options: 
        -h              Show help  
        -v              transigner version
        " 
        exit 1;
}

#########################################################

# Unset variables in case they are host assigned and thus inherited 
#unset -v reads transcripts out_dir aln_p em_p data_type mode

# Process options 
while getopts ":a:e:d:m:hv" option; do
    case $option in 
        a) a_flag=1
            aln_p="${OPTARG}"
            ;;
        e) e_flag=1
            em_p="${OPTARG}"
            ;;
        d) d_flag=1
            data_type="${OPTARG}"
            ;;
        m) m_flag=1
            mode="${OPTARG}"
            ;;
        h) Help
            ;;
        v) transigner --version >&2
            exit 1
            ;;
        \?) echo "Please enter in valid arguments, see -h for help" >&2 
            exit 1
            ;;
        :) echo "Option -${OPTARG} requires an argument" >&2
            exit 1
            ;;
    esac
done

# Need shift statement to handle positional args $1 $2 $3  
shift "$((OPTIND-1))"

#aln_p=24 # number of threads to use for minimap2 alignment
#em_p=4 # number of threads to use for EM iterations
#data_type="ont_drna" # input data type
#mode="default"

# Check if all mandatory options were provided
if [ $a_flag -eq 0 ] || [ $e_flag -eq 0 ] || [ $d_flag -eq 0 ] || [ $m_flag -eq 0 ] ; then
  echo "Error: All options -a, -e, -d and -m are mandatory."
  Help
fi

# Check if positional arguments are provided
if [ $# -lt 3 ]; then
  echo "Error: Missing positional arguments."
fi

# Now set variables for positional arguments 
reads="$1"
transcripts="$2"
out_dir="$3"

# Safety options
set -x

transigner align -q "$reads" -t "$transcripts" -d "$out_dir" \
   -o "aligned.bam" -p "$aln_p"
if [ "$mode" == "default" ]; then
    transigner pre -i "${out_dir}/aligned.bam" -d "$out_dir" 
    transigner em -s "${out_dir}/scores.csv" -d "$out_dir" \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype "$data_type" -p "$em_p"
elif [ "$mode" == "psw" ]; then
    transigner pre -i "${out_dir}/aligned.bam" -d "$out_dir" --use-psw
    transigner em -s "${out_dir}/scores.csv" -d "$out_dir" \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype "$data_type" -p "$em_p" --no-drop
elif [ "$mode" == "spiked" ]; then
    transigner pre -i "${out_dir}/aligned.bam" -d "$out_dir" --spiked
    transigner em -s "${out_dir}/scores.csv" -d "$out_dir" \
        -u "${out_dir}/unmapped.txt" -m "${out_dir}/tmap.csv" \
        -dtype "$data_type" -p "$em_p"
else
    echo "unrecognized mode"
fi

# END easy.sh 
