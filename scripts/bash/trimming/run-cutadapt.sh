#!/bin/bash

####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :



-i : input fastq file (= path to fastq read 1 if paired-end)

-o : output_path = path where you want to store the quantificaton results

-p : path to cutadapt pipeline

-l : sequencing library type, can be 'PE' (for paired end) or 'SE' (for single end) 

-args : other cutadapt parameters, e.g. the sequence of the adapter to trim, the minimum overlap : '-a ATTGGATCAATGCGGC -O 4'.'.

-emp : If set to 1, creates empty output files and folders. Not implemented yet, just exits when called.




"
        exit 0
fi



#######  ARGUMENTS #######

while [[ $# -gt 0 ]];
do
    opt="$1";
    shift; 
    case "$opt" in       
        -i) input_r1="$1"; shift;;            
        -o) output_path="$1"; shift;;
        -l) lib="$1"; shift;;
        -p) cutadapt="$1"; shift;;
        -args) args="$1"; shift;;
        -emp) emp="$1"; shift;;
    esac
done


####### IF A PATH TO CUTADAPT IS SET WITH -P, ADD IT TO $PATH #######

if [[ -n $cutadapt ]] || [[ $cutadapt != 'def' ]]
then
        export PATH="$cutadapt/bin:$PATH"
fi

### GET SAMPLE NAME

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}    



####### EMPTY MODE TO ONLY CREATE EMPTY OUTPUT FILES AND FOLDERS -- useful for qsub-all-fastq.sh #######

emp=${emp:-'0'}

if [[ $emp = 1 ]]
then
	if [[ $lib = PE ]]
        then
        	input_r2=${input_r1/1.fastq/2.fastq}
                touch $output_path/$sname/trimmed-$(basename $input_r1)
                touch $output_path/$sname/trimmed-$(basename $input_r2)
        elif [[ $lib = SE ]]
        then
                touch $output_path/$sname/trimmed-$(basename $input_r1)
 	fi
	exit 1
fi


####### CREATE LOG FILE FOR THIS SCRIPT #######

                                                                                    
mkdir -p $output_path/$sname/logs/ 
log_file=$output_path/$sname/logs/$(basename $0)_log.txt
touch $log_file

echo "$(date)

########################################################
############## Trim reads with cutadapt  ###############
########################################################

cutadapt version : $(cutadapt --version)


Run script    $(basename $0)    with args :

-i : $input_r1
-o : $output_path
-p : $cutadapt
-l : $lib
-args : ${args[@]}







" >> $log_file




####### RUN CUTADAPT #######

if [[ $lib = 'PE' ]]
        then			
		input_r2=${input_r1/1.fastq/2.fastq} 
                cutadapt ${args[@]} -o $output_path/$sname/trimmed-$(basename $input_r1) -p $output_path/$sname/trimmed-$(basename $input_r2) $input_r1 $input_r2
elif [[ $lib = 'SE' ]]
        then
                cutadapt ${args[@]} -o $output_path/$sname/trimmed-$(basename $input_r1) $input_r1
else
        echo "

        ERROR : seq parameter $lib is not valid. Can be either PE or SE." >> $log_file
fi

echo "











" >> $log_file
