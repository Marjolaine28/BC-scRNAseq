#!/bin/bash


if [[ $1 = '--help' ]]
        then
                echo "


--------------- HELP ---------------




Script    $(basename $0)    must be run with args :




-g : genome reference

-t : transcriptome reference

-o : output folder



"
exit 0
fi




while [[ $# -gt 0 ]];
do
    opt="$1";
    shift; 
    case "$opt" in
        -g) reference_genome="$1"; shift;;
        -t) reference_transcriptome="$1"; shift;;
        -o) output_folder="$1"; shift;;
    esac
done

# Create decoys and gentrome

mkdir -p $output_folder

grep "^>" <(zcat $reference_genome) | cut -d " " -f 1 > $output_folder/decoys.txt
sed -i -e 's/>//g' $output_folder/decoys.txt
cat $reference_transcriptome $reference_genome > $output_folder/gentrome.fa.gz
