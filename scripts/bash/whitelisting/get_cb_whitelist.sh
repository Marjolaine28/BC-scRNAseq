#!/bin/bash


####### HELP PANEL #######

if [[ $1 = '--help' ]]
then
    echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with a :



-f : input_r1 = path to fastq read 1 (that contains the barcodes)

-o : output_folder = path where you want to store the whitelist

-p : path to the python script that get top N barcodes and remove collisions among them (with connecte component procedure)

-args : other arguments to pass to the python script (e.g. '--top=3000')

-emp : If set to 1, creates empty output files and folders. Useful for qsub-all-fastq.sh when holding jobs.




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
        -p) python_whitelisting_script="$1"; shift;;
        -args) args="$1"; IFS=' ' declare -a 'args=($args)'; shift;;
        -emp) emp="$1"; shift;; 
    esac
done





### GET SAMPLE NAME

IFS='/' read -ra sname <<< $input_r1
sname=${sname[*]: -2:1}    




### DON'T OVERWRITE WHITELIST

if [[ -s $output_folder/$sname/whitelist.tsv ]]
then
    echo "Already performed whiteisting."
    exit 1
fi



####### EMPTY MODE TO ONLY CREATE EMPTY OUTPUT FILES AND FOLDERS -- useful for qsub-all-fastq.sh (holding jobs) #######

emp=${emp:-'0'}

if [[ $emp = 1 ]]
then
	touch $output_folder/$sname/whitelist.tsv
	exit 0
fi



####### ACTIVATE ENV #######

source ~/VirtualEns/python_3.9.5/bin/activate



####### WHITELISTING #######

if [[ ! -s $output_folder/$sname/raw_cb_freq.tsv ]]
then
    zcat $input_r1 | sed -n '2~4p'| awk '{print substr($0,0,12)}' | sort | uniq -c | sort -r -n > $output_folder/$sname/raw_cb_freq.tsv
    sed -i 's/^ *//' $output_folder/$sname/raw_cb_freq.tsv
    sed -i "s/ /\t/" $output_folder/$sname/raw_cb_freq.tsv
fi

python3 $python_whitelisting_script --input=$output_folder/$sname/raw_cb_freq.tsv --output=$output_folder/$sname/whitelist.tsv ${args[@]}






#####################  ADD LOG FILES ?? #######################
