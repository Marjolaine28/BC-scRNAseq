#!/bin/bash

####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with a :



-f : input_r1 = path to fastq read 1

-o : output_path = path where you want to store the quantificaton results

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
        -f) input_r1="$1"; shift;;            
        -o) output_path="$1"; shift;;
        -p) fastqc="$1"; shift;;
		-l) lib="$1"; shift;;
        -args) args="$1"; IFS=' ' declare -a 'args=($args)'; shift;;
        -emp) emp="$1"; shift;; 
    esac
done




### GET SAMPLE NAME

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}    



### DON'T OVERWRITE FASTQC REPORTS

if [[ -s $output_path/$sname/quality && ! -z "$(ls -A $output_path/$sname/quality)" ]]
then 
   echo "Already performed quality control of the reads."
   exit 1
fi



####### IF A PATH TO FASTQC IS SET WITH -P, ADD IT TO $PATH #######

if [[ -n $fastqc ]] || [[ $fastqc != 'def' ]]
then
    export PATH="$fastqc:$PATH"
fi



####### EMPTY MODE TO ONLY CREATE EMPTY OUTPUT FILES AND FOLDERS -- useful for qsub-all-fastq.sh #######
# No empty mode implemented ; not needed for now (no holding job waiting for it)

emp=${emp:-'0'}             # 0 by default

if [[ $emp = 1 ]]
then
	exit 0
fi



####### CREATE LOG FILE FOR THIS SCRIPT #######

                                                                                    
mkdir -p $output_path/${sname}/logs/ 
log_file=$output_path/${sname}/logs/$(basename $0)_log.txt
touch $log_file




echo "$(date)

###################################################################
############## Get fastQC sequencing quality report ###############
###################################################################

FastQC version : fastqc --v


Run script    $0    with args :

-f : $input_r1
-o : $output_path
-p : $fastqc
-l : $lib
-args : ${args[@]}







" >> $log_file




####### RUN FASTQC #######

mkdir -p $output_path/${sname}/quality

if [[ $lib = PE ]]
	then
		input_r2=${input_r1/1.fastq/2.fastq}
		fastqc $input_r1 --outdir=$output_path/${sname}/quality
		fastqc $input_r2 --outdir=$output_path/${sname}/quality
elif [[ $lib = SE ]]
	then
		fastqc $input_r1 --outdir=$output_path/${sname}/quality  
else
	echo "

	ERROR : seq parameter $lib is not valid. Can be either PE or SE." >> $log_file
fi

echo "











" >> $log_file

