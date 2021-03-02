#!/bin/bash

source ~/.bash_profile
module load torque

if [[ $1 = '--help' ]]
        then
                echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :




-w : walltime = time needed for the job, e.g. 12:00:00 for 12 hours

-m : mem = memory required for the job, e.g. 20GB

-d : decoys file

-t : transcriptome (or gentrome) fasta file

-k : k-mers length to build the index (can be a list) ; e.g. : '17 19 23 31'

-o : output folder 

-p : path to salmon pipeline

-args : other salmon args ; e.g. : '--gencode'



"
exit 0
fi


while [[ $# -gt 0 ]];
do
    opt="$1";
    shift; 
    case "$opt" in
        -w) walltime="$1"; shift;;
        -m) mem="$1"; shift;;
        -d) decoys="$1"; shift;;
        -t) transcriptome="$1"; shift;;
        -k) kmers="$1"; shift;;
        -o) output_folder="$1"; shift;;
        -p) salmon="$1"; shift;;
        -args) args="$1"; shift;;
    esac
done


# Default parameters 

walltime=${walltime:-'30:00:00'}
mem=${mem:-'100GB'}
kmers=${kmers:-"17 19 23 31"}

# Create index

IFS=' ' declare -a 'kmers=(${kmers[@]})'; 
for k in ${kmers[@]}
do
    mkdir -p $output_folder/index_k$k
    mkdir -p $output_folder/logs/
    if [[ -n $decoys ]]
    then
        echo "if [[ -n $salmon ]]; then export PATH="$PATH:$salmon/bin"; fi; salmon index -k $k -t $transcriptome -d $decoys -i $output_folder/index_k$k $args" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_folder/logs/ -N salmon_index_k${k}_decoys
    else
        echo "if [[ -n $salmon ]]; then export PATH="$PATH:$salmon/bin"; fi; salmon index -k $k -t $transcriptome -i $output_folder/index_k$k $args" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_folder/logs/ -N salmon_index_k${k}
    fi
done