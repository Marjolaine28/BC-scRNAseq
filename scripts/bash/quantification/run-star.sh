#!/bin/bash

module load star

####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with a :



-f : input_r1 = path to fastq read 1

-o : output_folder = path where you want to store the quantificaton results

-l : sequencing library type, can be 'PE' (for paired end), 'PEr2' (for paired end with cdna in read 2) or 'SE' (for single end) 

-p : path to star pipeline

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
        -o) output_folder="$1"; shift;;
        -l) lib="$1"; shift;;
        -p) star="$1"; shift;;
        -args) args="$1"; IFS=' ' declare -a 'args=($args)'; shift;;
        -emp) emp="$1"; shift;; 
    esac
done


####### EMPTY MODE TO ONLY CREATE EMPTY OUTPUT FILES AND FOLDERS -- useful for qsub-all-fastq.sh #######
# No empty mode implemented ; not needed for now (no holding job waiting for it)

emp=${emp:-'0'}             # 0 by default

if [[ $emp = 1 ]]
then
	exit 0
fi


####### IF A PATH TO STAR IS SET WITH -P, ADD IT TO $PATH #######

if [[ -n $star ]] || [[ $star != 'def' ]]
then
    export PATH="$star/bin:$PATH"
fi



####### GET SAMPLE NAME #######

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}


####### CHECK IF ANALYSIS WAS ALREADY PERFORMED #######

if [[ -f $output_folder/$sname/Aligned.sortedByCoord.out.bam ]]
then 
    echo "$output_folder/$sname/Aligned.sortedByCoord.out.bam already exists.
    
    
    "
    exit 1
fi
    

####### CREATE LOG FILE #######

mkdir -p $output_folder/$sname/logs/						                                                ### create log folder 
log_file=$output_folder/$sname/logs/$(basename $0)_log.txt
touch $log_file



echo "$(date)

############################################################
########### STAR quantification of RNA-seq reads ###########
############################################################

STAR version : $(STAR --v)


Run script    $(basename $0)    with args :

-f : $input_r1
-o : $output_folder
-p  $star
-args : ${args[@]}


(see also STAR log file in $output_folder/$sname/logs)" >> $log_file





####### RUN STAR #######

if [[ $lib = 'PE' ]]
    then
        input_r2=${input_r1/1.fastq/2.fastq}  	     			 ### get reads 2 files       
        STAR ${args[@]} --readFilesIn $input_r1 $input_r2 --readFilesCommand zcat --outFileNamePrefix $output_folder/$sname/
elif [[ $lib = 'SE' ]]
    then
        STAR ${args[@]} --readFilesIn $input_r1 --readFilesCommand zcat --outFileNamePrefix $output_folder/$sname/ 
elif [[ $lib = 'PEr2' ]]
    then
        input_r2=${input_r1/1.fastq/2.fastq}  	     			 ### get reads 2 files       
        STAR ${args[@]} --readFilesIn $input_r2 --readFilesCommand zcat --outFileNamePrefix $output_folder/$sname/
else
    echo "

    ERROR : lib parameter $lib is not valid. Can be either PE, PEr2 or SE." >> $log_file
fi

echo "











" >> $log_file
