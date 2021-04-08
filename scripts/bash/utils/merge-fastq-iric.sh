#!/bin/bash

path=$1 #path to raw data, e.g. $data/datasets/private/RNA-seq/bulk/DSP356/raw
lib=${2:-"PE"}

#for file in $path/Sample_*; do mv "$file" "${file#Sample_}";done;
#for file in $path/*; do mv "$file" "${file%_N*}";done;


samples=$(find $path -mindepth 1 -maxdepth 1 -type d)
for s in $samples
do

	echo "Merging $(basename $s) fastq files..."

	if [[ $lib == "PE" ]]
		then 
			n=$(find $s -maxdepth 1 -name "*1.fastq.gz" | wc -l)
			if [[ $n == 1 ]]
				then 
					mv $s/*1.fastq.gz $s/$(basename $s)_1.fastq.gz
					mv $s/*2.fastq.gz $s/$(basename $s)_2.fastq.gz
					echo "No need to merge.

					"
			else
					mkdir -p $s/unmerged
					mv $s/*fastq.gz $s/unmerged 2>/dev/null
					echo "$n files to merge ...

					"
					cat $(find $s/unmerged -name "*1.fastq.gz") > $s/$(basename $s)_1.fastq.gz
					cat $(find $s/unmerged -name "*2.fastq.gz") > $s/$(basename $s)_2.fastq.gz
			fi

	elif [[ $lib == "SE" ]]
		then 
			n=$(find $s -maxdepth 1 -name "*.fastq.gz" | wc -l)
			if [[ $n == 1 ]]
				then 
					mv $s/*.fastq.gz $s/$(basename $s).fastq.gz
					echo "No need to merge.

					"
			else
					mkdir -p $s/unmerged
					mv $s/*fastq.gz $s/unmerged 2>/dev/null
					echo "$n files to merge ...

					"
					cat $(find $s/unmerged -name "*.fastq.gz") > $s/$(basename $s).fastq.gz
			fi

	else
		echo "Uncorrect lib parameter."
	fi

done


