#!/bin/bash

input_path=$1 # path to downloaded fastqs
output_path=$2 # path to save the merged fastqs
depth=${3:-"6"} 
S=${5:-"all"} # samples to process
lib=${6:-"PE"} # PE for paired end sequencing, SE for single end




####### CACLCULATE TOTAL NUMBER OF SAMPLES (INCLUDING REPLICATES) THAT WILL BE PROCESSED #######

if [[ $S = 'all' ]]
then
    samples=($(find $input_path -mindepth $depth -maxdepth $(($depth+1)) -type d))
else
    samples=()
    for s in ${S[@]};
    do
        rep=($(find $input_path -mindepth $depth -maxdepth $(($depth+1)) -type d -name $s))
        samples=(${samples[@]} $rep)
    done
fi



for s in ${samples[@]}
do
 
	out=$output_path/$(basename $s)
	mkdir -p $out
	echo "Merging $(basename $s) fastq files..." 
	echo "Merging $(basename $s) fastq files..." >> $out/merge.log

	if [[ $lib == "PE" ]]
		then 
			if [[ -f $out/$(basename $s)_R1.fastq.gz ]]
			then
				echo "$(basename $s) already merged.
				
				" >> $out/merge.log
			else
				n1=($(find $s -name "*_R1*.fastq.gz"))
				n2=($(find $s -name "*_R2*.fastq.gz"))
				echo "${#n1[@]} files to merge ...
				
				"
				echo "${#n1[@]} files to merge ...

				" >> $out/merge.log
				if [[ ${#n1[@]} = 1 ]] ; then
					mv $n1 $out/$(basename $s)_R1.fastq.gz
					mv $n2 $out/$(basename $s)_R2.fastq.gz
				else
					cat ${n1[@]} > $out/$(basename $s)_R1.fastq.gz
					cat ${n2[@]} > $out/$(basename $s)_R2.fastq.gz
				fi
			fi

	elif [[ $lib == "SE" ]]
		then 
			if [[ -f $out/$(basename $s).fastq.gz ]]
			then
				echo "$(basename $s) already merged.
				
				" >> $out/merge.log
			else
				n=($(find $s -name "*.fastq.gz"))
				echo "${#n[@]} files to merge ...

				" >> $out/merge.log
				cat ${n[@]} > $out/$(basename $s).fastq.gz
			fi

	else
		echo "Uncorrect lib parameter."
	fi

done


