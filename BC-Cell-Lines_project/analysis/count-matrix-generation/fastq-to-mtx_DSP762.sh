#!/bin/bash

torque=${1:-"0"}
project_path=${2:-"$(git root)/BC-Cell-Lines_project"}
scripts_path=${3:-"$(git root)/scripts/bash"}


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

# wget --no-check-certificate -O - \
# "https://genomique.iric.ca/FastQList?key=2212-e9fcdc80cfd658683c786d89297a28fd&projectID=762&wget=1" \
# | wget --no-check-certificate -P $project_path/data/dsp762/downloaded -cri -

mkdir -p $project_path/data/dsp762/downloaded
cp -r /home/arion/davidm/TEST/dsp762/data/downloaded/* $project_path/data/dsp762/downloaded




############################## REORGANIZE AND RENAME FASTQ FILES ##############################

mkdir -p $project_path/data/dsp762/raw-fastqs

samples=$(find $project_path/data/dsp762/downloaded -mindepth 6 -maxdepth 7 -type d)
for s in $samples; do cp -r $s $project_path/data/dsp762/raw-fastqs; done
samples=$(find $project_path/data/dsp762/raw-fastqs -mindepth 1 -maxdepth 1 -type d)
for s in $samples; do d=$(basename $s); IFS='_' declare -a 'd=($d)'; mv $s $project_path/data/dsp762/raw-fastqs/${d[1]}; done           # this renaming step might be specific to the project (depends how the samples files were originally named by the IRIC genomic platform)




############################## MERGE FASTQ FILES ##############################

bash merge-fastq-iric.sh $project_path/data/dsp762/raw-fastqs




############################## QUALITY OF THE READS ##############################

bash $submit -w 03:00:00 -m 20gb -h "0" -r $scripts_path/quality/run-fastqc.sh -i $project_path/data/dsp762/raw-fastqs -o $project_path/data/dsp762/raw-fastqs -l "PE" -s "all" -p $project_path/pipelines/fastqc/FastQC-0.11.9




# ############################## TRIMMING READS ##############################

bash $submit -w 06:00:00 -m 30gb -h "0" -r $scripts_path/trimming/run-cutadapt.sh -i $project_path/data/dsp762/raw-fastqs -o $project_path/data/dsp762/trimmed-fastqs/cutadapt -l "PE" -s "all" -p $project_path/pipelines/cutadapt/cutadapt-3.2 -args '--times 8 -m 20:20 -a CTGTCTCTTATACACATCTC -A GTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGCCGTATC'




# ############################## QUALITY OF THE READS (TRIMMED) ##############################

bash $submit -w 03:00:00 -m 20gb -h "1" -r $scripts_path/quality/run-fastqc.sh -i $project_path/data/dsp762/trimmed-fastqs/cutadapt -o $project_path/data/dsp762/trimmed-fastqs/cutadapt -l "PE" -s "all" -p $project_path/pipelines/fastqc/FastQC-0.11.9




############################## SALMON INDEXING ##############################

output_folder_index=$project_path/pipelines/salmon/salmon-1.4.0/files/$ref/pc-decoys      # pc-decoys means that the index is build on the protein-coding gene set and using decoys

bash salmon-get-decoys.sh -g $project_path/references/$ref/genome/GRCh38.primary_assembly.genome.fa.gz \
-t $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz -o $output_folder_index

if [[ $torque = 0 ]]
then
    export PATH="$PATH:$project_path/pipelines/salmon/salmon-1.4.0/bin"
    salmon index -k 17 -t $output_folder_index/gentrome.fa.gz -d $output_folder_index/decoys.txt -i $output_folder_index/index_k17 --gencode
    # salmon index -k 31 -t $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz -i $output_folder_index/index_k31 --gencode

else
    bash qsub-salmon-indexing.sh -w 15:00:00 -m 150gb -t $output_folder_index/gentrome.fa.gz -d $output_folder_index/decoys.txt -o $output_folder_index -k 17 -p $project_path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
    # bash qsub-salmon-indexing.sh -w 3:00:00 -m 10gb -t $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz -o $output_folder_index -k 31 -p $project_path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
fi





############################## QUANTIFICATION ##############################

# Get transcripts to gene mappings
zcat $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/|/\t/g' $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv
sed -i 's/>//g' $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv

# Run Alevin using forceCells 7000 and noWhitelist (more straightforward than Alevin default specific \
# whitelisitng procedure, which often fails with initial knee estimation and is performing a final classification of the quality the cells which is confusing)
bash $submit -w 30:00:00 -m 100gb -h '1 2n' -r $scripts_path/quantification/run-alevin.sh -i $project_path/data/dsp762/trimmed-fastqs/cutadapt \
-o $project_path/data/dsp762/quant/alevin/$ref/trimmed-reads-cutadapt/pc-decoys/forceCells-7000-noWh/raw -seq "PE" -s "all" -p $project_path/pipelines/salmon/salmon-1.4.0 \
-args "-l ISR -i $output_folder_index/index_k17 --tgMap $project_path/references/$ref/transcriptome/gencode.v34.pc_transcripts_txp2gene.tsv \
--dropseq --dumpMtx --dumpFeatures --dumpUmiGraph --dumpCellEq --dumpBfh --dumpArborescences --forceCells 7000 --noWhitelist"