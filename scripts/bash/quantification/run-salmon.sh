#!/bin/bash


####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with a :



-f : input_r1 = path to fastq read 1

-o : output_folder = path where you want to store the quantificaton results

-l : sequencing library type, can be 'PE' (for paired end) or 'SE' (for single end) 

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
        -o) output_folder="$1"; shift;;
        -l) lib="$1"; shift;;
        -p) salmon="$1"; shift;;
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


####### IF A PATH TO SALMON IS SET WITH -P, ADD IT TO $PATH #######

if [[ -n $salmon ]] || [[ $salmon != 'def' ]]
then
    export PATH="$salmon/bin:$PATH"
fi


####### GET SAMPLE NAME #######

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}


####### CHECK IF ANALYSIS WAS ALREADY PERFORMED #######

if [[ -f $output_folder/$sname/quant.genes.sf ]]
then 
    echo "$output_folder/$sname/quant.genes.sf already exists.
    
    
    "
    exit 1
fi
    

####### CREATE LOG FILE #######

mkdir -p $output_folder/$sname/logs/						                                                ### create log folder 
log_file=$output_folder/$sname/logs/$(basename $0)_log.txt
touch $log_file


echo "$(date)

################################################################
########### Salmon quantification of RNA-seq reads ###########
################################################################

Salmon version : $(salmon --v)


Run script    $(basename $0)    with args :

-f : $input_r1
-o : $output_folder
-p  $salmon
-l : $lib
-args : ${args[@]}


(see also salmon log file in $output_folder/$sname/logs)" >> $log_file




####### RUN SALMON #######

if [[ $lib = 'PE' ]]
    then
        input_r2=${input_r1/1.fastq/2.fastq}  	     			 ### get reads 2 files      
        salmon quant ${args[@]} -1 $input_r1 -2 $input_r2 -o $output_folder/$sname/
elif [[ $lib = 'SE' ]]
    then
        salmon quant ${args[@]} -r $input_r1 -o $output_folder/$sname/
else
    echo "

    ERROR : lib parameter $lib is not valid. Can be either PE or SE." >> $log_file
fi


echo "











" >> $log_file

