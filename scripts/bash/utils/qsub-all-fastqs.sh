#!/bin/bash

module load torque

if [[ $1 = '--help' ]]
then
        echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :



-w : walltime = time needed for the job, e.g. 12:00:00 for 12 hours

-m : mem = memory required for the job, e.g. 20GB

-h : pids of the jobs that must complete before queueing this one (pids must be separted by a space ' ')

-r : script to submit with qsub, e.g. run-cutadapt.sh

-p : path to the pipeline called in the script, e.g. /home/arion/davidm/Sofwares/cutadapt/cutadapt-3.2

-f : input_path = path to the fastq file to process, e.g. /home/arion/davidm/Data/datasets/private/RNA-seq/sc/sc-MCF7_DSP779/raw

-o : output_path = path where you want to store output, e.g. /home/arion/davidm/Data/datasets/private/RNA-seq/sc/sc-MCF7_DSP779/trimmed

-l : lib = the sequencing type, can be 'PE' (for paired end) or 'SE' (for single end), or also 'ISR' (salmon nomenclature)

-s : samples = the samples / batch / condition / cell line that you want to analyse, e.g. 'ctrl er1' or 'GRHL2-treated GRHL2-non-treated' or 'all'

-args : additional arguments to pass to the script, e.g. adapter sequence for cutadapt or index for alevin.

-savepids : stores the pids in the path given (if savepids not passed as argument, pids are not stored)




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
        -savepids) savepids="$1"; shift;;
    esac
done



###### AUXILARY FUNCTIONS #######

arrayGet() { 
    local array=$1 index=$2
    local i="${array}_$index"
    printf '%s' "${!i}"
}

function join_by { local d=${1-} f=${2-}; if shift 2; then printf %s "$f" "${@/#/$d}"; fi; }



###### DEFAULT ARGS #######


lib=${lib:-'PE'}
hold=${hold:-''}
samples=${samples:-'all'}  
pipeline=${pipeline:-'def'}
savepids=${savepids:-'0'}



###### CONFIGURE TMP DIRCTORY #######

mkdir -p $savepids
if [[ -s $savepids/pids.txt ]]
then
    rm  $savepids/pids.txt;
fi


####### CACLCULATE TOTAL NUMBER OF SAMPLES (INCLUDING REPLICATES) THAT WILL BE PROCESSED #######

if [[ $samples = 'all' ]]
then
    if [[ $lib = SE ]]
    then
        S=($(find $input_path -maxdepth 2 -name "*.fastq*"))
    else
        S=($(find $input_path -maxdepth 2 -name "*R1.fastq*"))
    fi
else
    S=()
    for s in ${samples[@]};
    do
            if [[ $lib = SE ]]
            then
                    rep=($(find $input_path/$s/ -maxdepth 1 -name "*.fastq*"))
            else
                    rep=($(find $input_path/$s/ -maxdepth 1 -name "*R1.fastq*"))
                    S=(${S[@]} $rep)
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
-savepids : $savepids










" >> $log_file




####### SUBMIT JOBS #######


IFS=';'; hold=($hold);
run=$(basename $run_path);

echo ${hold[0]}
echo ${hold[1]}


for s in ${S[@]}                                                                                                            ### loop over the samples
do
    IFS='/' read -ra sname <<< $s
    sname=${sname[*]: -2:1}                                                                                                 ### get sample name
    mkdir -p $output_path/$sname/logs/						                                                                ### create log folder 
    $run_path -emp 1 -f $s -o $output_path -l $lib -p $pipeline; ec=$?                                                      ### it is useful to create empty files and folder for job with hold status (they take output of previous task as input, but if the output is not available yet the holding task would have no input)

    if [[ $ec == 1 ]]; 
    then
        echo "Job not submitted for $sname."                                                                                ### error code 1 is returned when the output already exists (avoid overwriting)
        continue
    fi

    pids="";
    IFS=';';                                                                                                     ### get the different jobs to wait for
    for h in ${hold[@]}; do
        echo $h
        IFS=','; holdx=($h);                                                                                                ### get the different sample-specific jobs to wait for
        if [[ ${#holdx[@]} > 1 ]];                                                                                          ### case only more than one sample
        then 
            IFS=' '; sx=(${holdx[0]}); pidx=(${holdx[1]});                                                                  ### get pid for the samples
            for index in "${!sx[@]}"; do
                if [[ ${sx[$index]} = $sname ]]
                then
                    pid=${pidx[$index]}
                    pids="$pids $pid";                                                                                      ### append the sample-specific job pid
                fi
            done
        else 
            pids="$pids $h";  
        fi; 
    done;
    IFS=' '; pids=($pids); pids=$(join_by : ${pids[@]});
    if [[ $pids = "" ]]                                                                                                     ### submit jobs ; hold option allows to wait for some previous task to finish before exection
    then
        pid_job=$(echo "$run_path -f $s -o $output_path -p $pipeline -args '$args' -l $lib" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_path/${sname}/logs/ -N "${run}_${sname}")
        echo "$pid_job submitted ($(basename $run_path) on $sname)."
    else
        echo "Wait jobs $pids to complete before queueing :"
        pid_job=$(echo "$run_path -f $s -o $output_path -p $pipeline -args '$args' -l $lib" | qsub -V -l nodes=1,mem=$mem,vmem=$mem,walltime=$walltime -j oe -d $output_path/${sname}/logs/ -N "${run}_${sname}" -W depend=afterok:$pids) 
        echo "$pid_job in hold status ($(basename $run_path) on $s)."

    fi

    if [[ ! $savepids = 0 ]]
    then
        mkdir -p $savepids
        echo $sname ${pid_job} $run_path >> $savepids/pids.txt
    fi
done