#!/bin/bash

source ~/.bash_profile

####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with a :



-i : input_r1 = path to fastq read 1

-o : output_folder = path where you want to store the quantificaton results

-p : path to salmon pipeline

-args : other alevin parameters (e.g. index, txp2gene file, --dropseq).

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
        -o) output_folder="$1"; shift;;
        -p) fastqc="$1"; shift;;
		-l) lib="$1"; shift;;
        -args) args="$1"; IFS=' ' declare -a 'args=($args)'; shift;;
        -emp) emp="$1"; shift;; 
    esac
done



####### EMPTY MODE TO ONLY CREATE EMPTY OUTPUT FILES AND FOLDERS -- useful for qsub-all-fastq.sh #######
# No empty mode implemented ; not needed for now (no holding job waiting for it)

emp=${emp:-'0'}             # 0 by default

if [[ $emp = 1 ]]
then
	exit 1
fi


####### IF A PATH TO FASTQC IS SET WITH -P, ADD IT TO $PATH #######

if [[ -n $fastqc ]] || [[ $fastqc != 'def' ]]
then
    export PATH="$PATH:$fastqc"
fi



####### CREATE LOG FILE FOR THIS SCRIPT #######

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}                                                                                     ### get sample name
mkdir -p $output_folder/$sname/logs/						                                                ### create log folder 
log_file=$output_folder/$sname/logs/$(basename $0)_log.txt
touch $log_file


echo "$(date)

###################################################################
############## Get fastQC sequencing quality report ###############
###################################################################

FastQC version : fastqc --v


Run script    $0    with args :

-i : $input_r1
-o : $output_folder
-p : $fastqc
-l : $lib
-args : ${args[@]}







" >> $log_file




####### RUN FASTQC #######

if [[ $lib = PE ]]
	then
		input_r2=${input_r1/1.fastq/2.fastq}
		fastqc $input_r1 --outdir=$output_folder/${sname}/quality
		fastqc $input_r2 --outdir=$output_folder/${sname}/quality
elif [[ $lib = SE ]]
	then
		fastqc $input_r1 --outdir=$output_folder/${sname}/quality  
else
	echo "

	ERROR : seq parameter $lib is not valid. Can be either PE or SE." >> $log_file
fi

echo "











" >> $log_file

