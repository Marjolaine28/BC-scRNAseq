#!/bin/bash

torque=${1:-"0"}
data_path=${2:-"$(git root)/BC-Cell-Lines_project/data/dsp992"}
scripts_path=${3:-"$(git root)/scripts/bash"}
references_path=${3:-"$(git root)/references"}
pipelines_path=${3:-"$(git root)/pipelines"}




############################## PREPARE WORKING ENVIRONNEMENT ##############################

source ~/.bash_profile

export PATH="$PATH:$scripts_path/utils/" 
export PATH="$PATH:$scripts_path/quantification/" 

if [[ $torque = 0 ]]
then
    submit=run-all-fastqs.sh
else
    submit=qsub-all-fastqs.sh
fi

ref=human/assembly__GRCh38-hg38/annotation__gencode/gencode_34




############################## DOWNLOAD FASTQ FILES ##############################

wget --no-check-certificate -O - \
"https://genomique.iric.ca/FastQList?key=2212-e9fcdc80cfd658683c786d89297a28fd&projectID=992&wget=1" \
| wget --no-check-certificate -P $data_path/downloaded -cri -




############################## REORGANIZE AND RENAME FASTQ FILES ##############################

mkdir -p $data_path/raw-fastqs

samples=$(find $data_path/downloaded -mindepth 6 -maxdepth 7 -type d)
for s in $samples; do cp -r $s $data_path/raw-fastqs; done
samples=$(find $data_path/raw-fastqs -mindepth 1 -maxdepth 1 -type d)
for s in $samples; do d=$(basename $s); IFS='_' declare -a 'd=($d)'; mv $s $data_path/raw-fastqs/${d[0]}; done           # this renaming step might be specific to the project (depends how the samples files were originally named by the IRIC genomic platform)




############################## MERGE FASTQ FILES ##############################

merge-fastq-iric.sh $data_path/raw-fastqs




############################## QUALITY OF THE READS ##############################

$submit -w 03:00:00 -m 20gb -h "0" -r $scripts_path/quality/run-fastqc.sh -i $data_path/raw-fastqs -o $data_path/raw-fastqs -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9




# ############################## TRIMMING READS ##############################

$submit -w 06:00:00 -m 30gb -h "0" -r $scripts_path/trimming/run-cutadapt.sh -i $data_path/raw-fastqs -o $data_path/trimmed-fastqs/cutadapt -l "PE" -s "all" -p $pipelines_path/cutadapt/cutadapt-3.2 -args '--times 8 -m 20:20 -a CTGTCTCTTATACACATCTC -A GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATC'




# ############################## QUALITY OF THE READS (TRIMMED) ##############################

$submit -w 03:00:00 -m 20gb -h "1" -r $scripts_path/quality/run-fastqc.sh -i $data_path/trimmed-fastqs/cutadapt -o $data_path/trimmed-fastqs/cutadapt -l "PE" -s "all" -p $pipelines_path/fastqc/FastQC-0.11.9




############################## SALMON INDEXING ##############################

output_folder_index=$pipelines_path/salmon/salmon-1.4.0/files/$ref/pc-decoys      # pc-decoys means that the index is build on the protein-coding gene set and using decoys

salmon-get-decoys.sh -g $references_path/$ref/genome/GRCh38.primary_assembly.genome.fa.gz \
-t $references_path/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz -o $output_folder_index

if [[ $torque = 0 ]]
then
    export PATH="$pipelines_path/salmon/salmon-1.4.0/bin:$PATH"
    salmon index -k 17 -t $output_folder_index/gentrome.fa.gz -d $output_folder_index/decoys.txt -i $output_folder_index/index_k17 --gencode

else
    qsub-salmon-indexing.sh -w 15:00:00 -m 150gb -t $output_folder_index/gentrome.fa.gz -d $output_folder_index/decoys.txt -o $output_folder_index -k 17 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
fi





############################## QUANTIFICATION ##############################

# Get transcripts to gene mappings
zcat $references_path/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $references_path/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $references_path/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/>//g' $references_path/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv

# Run Alevin using forceCells 7000 and noWhitelist (more straightforward than Alevin default specific \
# whitelisitng procedure, which often fails with initial knee estimation and is performing a final classification of the quality the cells which is confusing)
$submit -w 30:00:00 -m 100gb -h '1 2n' -r $scripts_path/quantification/run-alevin.sh -i $data_path/trimmed-fastqs/cutadapt \
-o $data_path/quant/alevin/$ref/trimmed-reads-cutadapt/pc-decoys-k17-1.4.0/forceCells-7000-noWh/raw -seq "PE" -s "all" -p $pipelines_path/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/index_k17 --tgMap $references_path/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv \
--dropseq --dumpMtx --dumpFeatures --dumpUmiGraph --dumpCellEq --dumpBfh --dumpArborescences --forceCells 7000 --noWhitelist"