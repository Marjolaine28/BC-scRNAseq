#!/bin/bash

input_path=$1 #path to raw data, e.g. $data/datasets/private/RNA-seq/bulk/DSP356/raw
output_path=$2
lib=${3:-"PE"}

#for file in $path/Sample_*; do mv "$file" "${file#Sample_}";done;
#for file in $path/*; do mv "$file" "${file%_N*}";done;

samples=$(find $input_path -mindepth 6 -maxdepth 7 -type d)
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
				n1=($(find $s -name "*_R1_*.fastq.gz"))
				n2=($(find $s -name "*_R2_*.fastq.gz"))
				echo "${#n1[@]} files to merge ...

				" >> $out/merge.log
				cat ${n1[@]} > $out/$(basename $s)_R1.fastq.gz
				cat ${n2[@]} > $out/$(basename $s)_R2.fastq.gz
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


