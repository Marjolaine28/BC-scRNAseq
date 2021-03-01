#!/bin/bash

project_path=${1:-"../.."}
scripts_path=${2:-"../../../scripts/"}

source ~/.bash_profile

mkdir -p $project_path/data
mkdir -p $project_path/pipelines
mkdir -p $project_path/references


############################## INSTALL REQUIRED PIPELINES ##############################

## FastQC 0.11.8
wget https://github.com/s-andrews/FastQC/archive/v0.11.8.tar.gz -P $project_path/pipelines/fastqc
tar -xvf $project_path/pipelines/fastqc/v0.11.8.tar.gz -C $project_path/pipelines/fastqc
$project_path/pipelines/fastqc/v0.11.8.tar.gz 

## Cutadapt 3.2
mkdir -p $project_path/pipelines/cutadapt
python3 -m venv $project_path/pipelines/cutadapt/cutadapt-3.2
$project_path/pipelines/cutadapt/cutadapt-3.2/bin/pip install --upgrade pip
$project_path/pipelines/cutadapt/cutadapt-3.2/bin/pip install cutadapt==3.2
# export PATH="$PATH:$project_path/pipelines/cutadapt/cutadapt-3.2/bin"

## Salmon 1.4.0
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz -P $project_path/pipelines/salmon 
tar -xvf $project_path/pipelines/salmon/salmon-1.4.0_linux_x86_64.tar.gz -C $project_path/pipelines/salmon
rm $project_path/pipelines/salmon/salmon-1.4.0_linux_x86_64.tar.gz
mv $project_path/pipelines/salmon/salmon-latest_linux_x86_64/ $project_path/pipelines/salmon/salmon-1.4.0
# export PATH="$PATH:$project_path/pipelines/salmon/salmon-1.4.0/bin"




############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################

ref=human/assembly__GRCh38-hg38/annotation__gencode/gencode_34
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_transcripts.fa.gz -P $project_path/references/$ref/transcriptome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $project_path/references/$ref/genome



############################## DOWNLOAD FASTQ FILES ##############################

wget --no-check-certificate -O - \
"https://genomique.iric.ca/FastQList?key=2212-e9fcdc80cfd658683c786d89297a28fd&projectID=762&wget=1" \
| wget --no-check-certificate -P $project_path/data/downloaded -cri -


############################## REORDER AND RENAME FASTQ FILES ##############################

mkdir -p $project_path/data/raw-fastqs
for dir in $project_path/data/downloaded/*/*/*/*/*/*; do cp -r $dir $project_path/data/raw-fastqs; done
for dir in $project_path/data/raw-fastqs/*; do IFS='_' declare -a 's=($dir)'; mv $dir ${s[@]: -1:2}; done           # this renaming step might be specific to the project (depends how the samples files were originally named by the IRIC genomic platform)


############################## MERGE FASTQ FILES ##############################

bash merge-fastq-iric.sh $project_path/raw-fastqs/


############################## QUALITY OF THE READS ##############################

bash qsub-all-fastq.sh -w 03:00:00 -m 20gb -h "1" -r run-fastqc.sh -i $project_path/raw-fastqs -o $project_path/raw-fastqs -l "PE" -s "all" -p $project_path/pipelines/fastqc/FastQC-0.11.8


############################## TRIMMING READS ##############################

bash qsub-all-fastq.sh -w 06:00:00 -m 30gb -h "0" -r run-cutadapt.sh -i $project_path/raw-fastqs -o $project_path/trimmed/cutadapt -l "PE" -s "all" -p $project_path/pipelines/cutadapt/cutadapt-3.2 -args '--times 8 -m 20:20 -a CTGTCTCTTATACACATCTC -A GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATC'



############################## QUALITY OF THE READS (TRIMMED) ##############################

bash qsub-all-fastq.sh -w 03:00:00 -m 20gb -h "1" -r run-fastqc.sh -i $project_path/trimmed/cutadapt -o $project_path/trimmed/cutadapt -l "PE" -s "all" -p $project_path/pipelines/fastqc/FastQC-0.11.8


############################## SALMON INDEXING ##############################

output_folder_index=$project_path/pipelines/salmon/salmon-1.4.0/files/$ref/pc-decoys      # pc-decoys means that the index is build on the protein-coding gene set and using decoys
bash salmon-get-decoys.sh -g $project_path/references/$ref/genome/GRCh38.primary_assembly.genome.fa.gz \
-t $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz -o $output_folder_index
bash qsub-salmon-indexing.sh -w 15:00:00 -m 150gb -t $output_folder_index/gentrome.fa.gz -d $output_folder_index/decoys.txt -o $output_folder_index -k 17 -p $project_path/pipelines/salmon/salmon-1.4.0 -args '--gencode'



############################## QUANTIFICATION ##############################

# Get transcripts to gene mappings
zcat $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/>//g' $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv

# Run Alevin using forceCells 7000 and noWhitelist (more straightforward than Alevin default specific \
# whitelisitng procedure, which often fails with initial knee estimation and is performing a final classification of the quality the cells which is confusing)
bash qsub-all-fastq.sh -w 30:00:00 -m 100gb -h '1 2n' -r run-alevin.sh -i $project_path/trimmed/cutadapt \
-o $project_path/quant/alevin/$ref/trimmed-cutadapt/pc-decoys/forceCells-7000-noWh -seq "PE" -s "all" -p $project_path/pipelines/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/index_k17 --tgMap $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv \
--dropseq --dumpMtx --dumpFeatures --dumpUmiGraph --dumpCellEq --dumpBfh --dumpArborescences --forceCells 7000 --noWhitelist"