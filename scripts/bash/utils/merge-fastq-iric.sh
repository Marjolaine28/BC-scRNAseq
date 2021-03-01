#!/bin/bash

path=$1 #path to raw data, e.g. $data/datasets/private/RNA-seq/bulk/DSP356/raw

#for file in $path/Sample_*; do mv "$file" "${file#Sample_}";done;
#for file in $path/*; do mv "$file" "${file%_N*}";done;

for s in $path/*
do
	mkdir -p $s/unmerged
	mv $s/*fastq.gz $s/unmerged 2>/dev/null	
	cat $(find $s/unmerged -name "$(basename $s)*_R1_*") > $s/$(basename $s)_1.fastq.gz
	cat $(find $s/unmerged -name "$(basename $s)*_R2_*") > $s/$(basename $s)_2.fastq.gz
	done
done

