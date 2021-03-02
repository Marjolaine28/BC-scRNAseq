#!/bin/bash

if [[ $1 = '--help' ]]
then
        echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :



-r : script to run, e.g. run-cutadapt.sh

-p : path to the pipeline called in the script, e.g. /home/arion/davidm/Sofwares/cutadapt/cutadapt-3.2

-i : input_path = path to the data to process, e.g. /home/arion/davidm/Data/datasets/private/RNA-seq/sc/sc-MCF7_DSP779/raw

-o : output_path = path where you want to store output, e.g. /home/arion/davidm/Data/datasets/private/RNA-seq/sc/sc-MCF7_DSP779/trimmed

-l : lib = the sequencing type, can be 'PE' (for paired end) or 'SE' (for single end), or also 'ISR' (salmon nomenclature)

-s : samples = the samples / batch / condition / cell line that you want to analyse, e.g. 'ctrl er1' or 'GRHL2-treated GRHL2-non-treated' or 'all'

-args : additional arguments to pass to the script, e.g. adapter sequence for cutadapt or index for alevin.




"
exit 0
fi


#######  ARGUMENTS #######

while [[ $# -gt 0 ]];
do
    opt="$1";
    shift; 
    case "$opt" in
	    -r) run_path="$1"; shift;;  
        -p) pipeline="$1"; shift;;         
        -i) input_path="$1"; shift;;            
        -o) output_path="$1"; shift;;
	    -l) lib="$1"; shift;;
        -s) samples="$1"; IFS=' ' declare -a 'samples=($samples)'; shift;;
        -args) args="$1"; shift;;
    esac
done

lib=${lib:-'PE'}
samples=${samples:-'all'}  
pipeline=${pipeline:-'def'} 


####### CACLCULATE TOTAL NUMBER OF SAMPLES (INCLUDING REPLICATES) THAT WILL BE PROCESSED #######

if [[ $samples = 'all' ]]
then
    if [[ $lib = PE ]]
    then
        S=($(find $input_path -maxdepth 2 -name *1.fastq*))
    elif [[ $lib = SE ]]
    then
        S=($(find $input_path -maxdepth 2 -name *.fastq*))
    fi
else
    S=()
    for s in ${samples[@]};
    do
            if [[ $lib = PE ]]
            then
                    rep=($(find $input_path/$s/ -maxdepth 1 -name *1.fastq*))
                    S=(${S[@]} $rep)
            elif [[ $lib = SE ]]
            then
                    rep=($(find $input_path/$s/ -maxdepth 1 -name *.fastq*))
            fi
    done
fi
tot=${#S[@]}



###### CREATE A LOG FILE FOR THIS SCRIPT #######

mkdir -p $output_path/logs

log_file=$output_path/logs/$(basename $0)_log.txt
touch $log_file

echo "$(date)

$tot jobs submited.

Run script    $(basename $0)    with args :

-p : $pipeline
-r : $run_path
-i : $input_path
-o : $output_path
-l : $lib
-s : $samples
-args : ${args:-"default"}










" >> $log_file




####### SUBMIT JOBS #######

run=$(basename $run_path) 

for s in ${S[@]}
do
        IFS='/' read -ra sname <<< $s
        sname=${sname[*]: -2:1}                                                                                                 ### get sample name
        mkdir -p $output_path/$sname/logs/						                                                                ### create log folder 
        # find  $output_path/$sname/ -mindepth 1 -maxdepth 1 -name '*fastq*' -exec rm -r 2>/dev/null "{}" \;	                ### remove existing results
        $run_path -i $s -o $output_path -p $pipeline -args "$args" -l $lib
done