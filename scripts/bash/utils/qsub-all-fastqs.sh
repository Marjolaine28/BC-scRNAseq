#!/bin/bash

module load torque

if [[ $1 = '--help' ]]
then
        echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :



-w : walltime = time needed for the job, e.g. 12:00:00 for 12 hours

-m : mem = memory required for the job, e.g. 20GB

-h : n when the job must wait for the n-th previous task to complete (hold status), else 0

-r : script to submit with qsub, e.g. run-cutadapt.sh

-p : path to the pipeline called in the script, e.g. /home/arion/davidm/Sofwares/cutadapt/cutadapt-3.2

-f : input_path = path to the fastq file to process, e.g. /home/arion/davidm/Data/datasets/private/RNA-seq/sc/sc-MCF7_DSP779/raw

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
        -w) walltime="$1"; shift;;        
        -m) mem="$1"; shift;;            
        -h) hold="$1"; shift;;
	-r) run_path="$1"; shift;;  
        -p) pipeline="$1"; shift;;         
        -f) input_path="$1"; shift;;            
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

-w : $walltime
-m: $mem
-h : $hold
-r : $run_path
-p : $pipeline
-f : $input_path
-o : $output_path
-l : $lib
-s : $samples
-args : ${args:-"default"}










" >> $log_file




####### SUBMIT JOBS #######

run=$(basename $run_path) 


if [[ $hold != 0 ]]
then
        IFS=' ' declare -a 'hold=($hold)'
        J=''
        for hi in ${hold[@]}
        do 
                if [[ ${hi:1:2} != n ]] && [[ ${hi:0:1} != 0 ]]
                then    
                        n=${hi:0:1}
                        c=0
                        for hj in ${hold[@]}
                        do
                                if [[ ${hj:1:2} = n ]] && [[ ${hj:0:1} != 0 ]] && [[ ${hj:0:1} < ${hi:0:1} ]]
                                then
                                        n=$(($n-1))
                                        c=$(($c+$tot))
                                fi
                        done
                        J=$J:$(qstat -u $USER | tail -n $(($n+$c)) | head -n 1 | awk '{print $1}' | cut -d"." -f1)
                fi
        done
fi


for s in ${S[@]}
do
        IFS='/' read -ra sname <<< $s
        sname=${sname[*]: -2:1}                                                                                                 ### get sample name
        mkdir -p $output_path/$sname/logs/						                                        ### create log folder 
        # find  $output_path/$sname/ -mindepth 1 -maxdepth 1 -name '*fastq*' -exec rm -r 2>/dev/null "{}" \;	                ### remove existing results
        $run_path -emp 1 -f $s -o $output_path -l $lib -p $pipeline                                                        ### it is useful to create empty files and folder for job with hold status (they take output of previous task as input, but if the output is not available yet the holding task would have no input)
		
        if [[ $hold = 0 ]]                                                                                                      ### submit jobs ; hold option allows to wait for some previous task to finish before exection
        then    
                echo "$run_path -f $s -o $output_path -p $pipeline -args '$args' -l $lib" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_path/${sname}/logs/ -N "${run}_${sname}"
        else	
                W=$J				
                for hj in ${hold[@]}
                do 
                        if [[ ${hj:1:2} = n ]] && [[ ${hj:0:1} != 0 ]]
                        then
                                n=${hj:0:1}
                                c=0
                                for hi in ${hold[@]}
                                do
                                        if [[ ${hi:1:2} != n ]] && [[ ${hi:0:1} != 0 ]] && [[ ${hi:0:1} < ${hj:0:1} ]]
                                        then
                                                n=$(($n-1))
                                                c=$(($c+1))
                                        fi
                                done
                                W=$W:$(qstat -u $USER | tail -n $(($tot*$n+$c)) | head -n 1 | awk '{print $1}' | cut -d"." -f1)
                        fi
                done
                W=${W:1}
                echo "Wait job $W to complete before queueing :"
                echo "$run_path -f $s -o $output_path -p $pipeline -args '$args' -l $lib" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_path/${sname}/logs/ -N "${run}_${sname}" -W depend=afterok:$W
                echo "(Hold status)"
        fi
done