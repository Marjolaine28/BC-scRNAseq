#!/bin/bash

torque=${1:-"0"}
data_path=${2:-"$(git root)/data"}
scripts_path=${3:-"$(git root)/scripts"}
pipelines_path=${4:-"$(git root)/pipelines"}
key=${5:-"2212-e9fcdc80cfd658683c786d89297a28fd"}







#################################################################################
######################### PREPARE WORKING ENVIRONNEMENT #########################
#################################################################################


export PATH="$scripts_path/bash/utils/:$PATH" 
export PATH="$scripts_path/bash/quantification/:$PATH"
export PATH="$pipelines_path/sratoolkit/sratoolkit.3.0.0-centos_linux64/bin/:$PATH"
export PATH="$pipelines_path/salmon/salmon-1.4.0/bin/:$PATH"
export PATH="$pipelines_path/star/2.7.10a/bin/:$PATH"
export PATH="$pipelines_path/samtools/samtools-1.16/bin/:$PATH"
export PATH="$pipelines_path/rsem/RSEM-1.3.3/:$PATH"



if [[ $torque = 0 ]]
then
    submit=run-all-fastqs.sh
else
    module load torque
    submit=qsub-all-fastqs.sh
fi


assembly=human/assembly__GRCh38-hg38

arrayGet() { 
    local array=$1 index=$2
    local i="${array}_$index"
    printf '%s' "${!i}"
}








#############################################################################
############################## DOWNLOAD FASTQS ##############################  
#############################################################################


# IRIC projets #
################

project_IDs=(356 280)


for p in ${project_IDs[@]}; do
    wget --no-check-certificate -O - \
    "https://genomique.iric.ca/FastQList?key=${key}&projectID=$p&wget=1" \
    | wget --no-check-certificate -P $data_path/iric/sc/dsp$p/downloaded -cri -

    merge-fastq-iric.sh $data_path/iric/sc/dsp$p/downloaded $data_path/iric/sc/dsp$p/raw-fastqs
done



# 56 cell lines #
#################

mkdir $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded

x=$(cat $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SRR_Acc_List.txt) ; 
for d in $x; do 
    prefetch $d ; 
    fasterq-dump --split-files $d; 
    rm -r $d; 
    gzip ${d}_*.fastq; 
done

sed -i -e 's/T47D Kbluc/T47D_Kbluc/g' $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/MDAMB175/MDAMB175VII/g' $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/SUM1315/SUM1315MO2/g' $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt

x=($(echo $(cut -f1 -d',' $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))
y=($(echo $(cut -f8 -d',' $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))

for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then mkdir -p ${y[$i]} ; mv ${x[$i]}*1.fastq.gz ${y[$i]}/${x[$i]}*R1.fastq.gz ; mv ${x[$i]}*2.fastq.gz ${y[$i]}/${x[$i]}*R2.fastq.gz ; fi ; done
for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then echo ${y[$i]} ; echo ${x[$i]} ; fi ; done


merge-fastq-iric.sh $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/raw-fastqs








#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################


if [[ ! -f $data_path/references/$assembly/GRCh38.primary_assembly.genome.fa.gz ]]; 
then 
    mkdir -p $data_path/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $data_path/references/$assembly
else
    echo "Reference already $assembly available."
fi


annot=annotation__gencode/gencode_34

if [[ ! -f $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz ]]; 
then 
    mkdir -p $data_path/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $data_path/references/$assembly/$annot
else
    echo "Reference already $annot available."
fi


annot=annotation__gencode/gencode_37

if [[ ! -f $data_path/references/$assembly/$annot/gencode.v37.transcripts.fa.gz ]]; 
then 
    mkdir -p $data_path/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v37.transcripts.fa.gz -P $data_path/references/$assembly/$annot
else
    echo "Reference already $annot available."
fi





##################################################################################
######################## GET TRANSCRIPTS TO GENES MAPPING ########################
##################################################################################


annot=annotation__gencode/gencode_34

if [[ ! -f $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv ]]; then
    zcat $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/|/\t/g' $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/>//g' $data_path/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
fi


annot=annotation__gencode/gencode_37

if [[ ! -f $data_path/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv ]]; then
    zcat $data_path/references/$assembly/$annot/gencode.v37.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $data_path/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
    sed -i 's/|/\t/g' $data_path/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
    sed -i 's/>//g' $data_path/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
fi






#############################################################################
############################## SALMON INDEXING ##############################    
#############################################################################


annot=annotation__gencode/gencode_34
output_folder_index=$pipelines_path/salmon/salmon-1.4.0/files/$assembly/$annot 



# Gets decoys with gencode 34 #
###############################

mkdir -p $output_folder_index/decoys

if [[ ! -f $output_folder_index/decoys/decoys.txt && ! -f $output_folder_index/decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $data_path/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys
fi 



# Indexing with decoys & gencode 34 #
#####################################

if [[ -d $output_folder_index/decoys/index_k31 ]]; then
    echo "decoys/index_k31 already generated.
    
    "
else
    if [[ $torque = 0 ]]; then
        export PATH="$pipelines_path/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 31 -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -i $output_folder_index/decoys/index_k31 --gencode

    else
        qsub-salmon-indexing.sh -w 5:00:00 -m 150gb -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -o $output_folder_index/decoys -k 31 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_decoys_k31=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)

    fi
fi




# Indexing without decoys & gencode 34 #
########################################

mkdir -p $output_folder_index/no-decoys

if [[ -d $output_folder_index/no-decoys/index_k31 ]]; then
    echo "no-decoys/index_k31 already generated.
    
    "
else
    if [[ $torque = 0 ]]
    then
        salmon index -k 31 -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -i $output_folder_index/no-decoys/index_k31 --gencode
    else
        qsub-salmon-indexing.sh -w 1:00:00 -m 10gb -t $data_path/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys -k 31 -p $pipelines_path/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_no_decoys_k31=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi







###########################################################################
############################## STAR INDEXING ##############################    
###########################################################################


# Create STAR index with GRCh38-hg38 assembly and gencode 37 annotation #
#########################################################################

annot=annotation__gencode/gencode_37
output_folder_index=$pipelines_path/star/2.7.10a/files/$assembly/$annot


echo "STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $output_folder_index/index \
--genomeFastaFiles $data_path/references/$assembly/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile $data_path/references/$assembly/$annot/gencode.v37.annotation.gtf" \
| qsub -V -l nodes=1,mem=200gb,vmem=200gb,walltime=48:00:00 -j oe -d $pipelines_path/star/2.7.10a/files/$assembly/$annot/index/logs/ -N "STAR-index-hg38-v37"









###########################################################################
############################## RSEM INDEXING ##############################    
###########################################################################


# Create RSEM reference with Gencode 37
#######################################

annot=annotation__gencode/gencode_37
output_folder_index=$pipelines_path/rsem/RSEM-1.3.3/files/$assembly/$annot 

rsem-prepare-reference --gtf $data_path/references/$assembly/$annot/gencode.v37.annotation.gtf \
$data_path/references/$assembly/GRCh38.primary_assembly.genome.fa $output_folder_index/index/ ####### QSUB ? WAIT PID ??









#############################################################################################################
##############################  MAPPING BULK DATA WITH SALMON FOR CORRELATIONS ##############################
#############################################################################################################


# 56 cells lines - decoys k31 & gencode 34 #
############################################

annot=annotation__gencode/gencode_34
output_folder_index=$pipelines_path/salmon/salmon-1.4.0/files/$assembly/$annot/index_k31 

qsub-all-fastqs.sh -w 20:00:00 -m 100gb -h 0 -r $scripts_path/bash/quantification/run-salmon.sh -f $data_path/public/bulk/BC-Cell-Lines-Panel_GSE48213/raw-fastqs \
-o $data_path/public/RNA-seq/bulk/BC-Cell-Lines-Panel_GSE48213/quant/salmon/$assembly/$annot/raw-reads/pc-decoys-k31-1.4.0/raw-counts -l "PE" -s "all" -p $pipelines_path/RNA-seq/salmon/salmon-1.4.0 \
-args "-l A -i $output_folder_index -g $data_path/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"


# MCF7 (DSP356) - no decoys k31 & gencode 34 #
##############################################


qsub-all-fastqs.sh -w 10:00:00 -m 60gb -r $scripts_path/bash/quantification/run-salmon.sh -f $data_path/iric/bulk/dsp356/raw-fastqs \
-o $data_path/iric/bulk/dsp356/quant/salmon/$assembly/$annot/raw-reads/no-decoys-k31-1.4.0/raw-counts -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
-p $pipelines_path/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/no-decoys/index_k31 \
-g $output_folder_index -g $data_path/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"


# MCF7 (DSP356) - decoys k31 & gencode 34 #
###########################################


qsub-all-fastqs.sh -w 10:00:00 -m 60gb -r $scripts_path/bash/quantification/run-salmon.sh -f $data_path/iric/bulk/dsp356/raw-fastqs \
-o $data_path/iric/bulk/dsp356/quant/salmon/$assembly/$annot/raw-reads/decoys-k31-1.4.0/raw-counts -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
-p $pipelines_path/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/decoys/index_k31 \
-g $output_folder_index -g $data_path/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"








####################################################################################################
############################## STAR ALIGNMENTS FOR CORRELATIONS & IGV ##############################
####################################################################################################


annot=annotation__gencode/gencode_37


# Align sc T47D (DSP762) in bulk mode with STAR (just for IGV, so no need to quantify after)
############################################################################################


qsub-all-fastqs.sh -w 30:00:00 -m 100gb -r $scripts_path/bash/quantification/run-star.sh -f $data_path/iric/sc/dsp762/trimmed-fastqs \ #### CHECK TRIMMED PATH
-o $data_path/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt -l "PEr2" -s "Sample_T47D" \
-args "--genomeDir $pipelines_path/star/2.7.10a/files/$assembly/$annot/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"





# Align bulk T47D (DSP280) with STAR
####################################


qsub-all-fastqs.sh -w 30:00:00 -m 100gb -h 0 -r $scripts_path/bash/quantification/run-star.sh -f $data_path/iric/bulk/dsp280/raw-fastqs \
-o $data_path/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads -l "PE" -s "all" \
-args "--genomeDir $pipelines_path/star/2.7.10a/files/$assembly/$annot/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"




### + Alignments already launched by the bioinfo plateforme  :

project_IDs=(550 1111)

for p in ${project_IDs[@]}; do
    wget --no-check-certificate -O - \
    "https://genomique.iric.ca/BamList?key=${key}&projectID=$p&wget=1" \
    | wget --no-check-certificate -P $data_path/iric/bulk/dsp$p/downloaded-align/star/$assembly/$annot/raw-reads -cri -
done

mv $data_path/iric/bulk/dsp$p/downloaded-align/star/$assembly/$annot/raw-reads/bioinfo*/*/*/*/* $data_path/iric/bulk/dsp$p/downloaded-align/star/$assembly/$annot/raw-reads
rm -r $data_path/iric/bulk/dsp$p/downloaded-align/star/$assembly/$annot/raw-reads/bioinfo*







####################################################################################################
############################## RSEM QUANTIFICATION FOR CORRELATIONS & IGV ##########################
####################################################################################################


annot=annotation__gencode/gencode_37


# RSEM quant for bulk T47D (DSP280)
###################################

x=($(echo $(find $data_path/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads -mindepth 1 -maxdepth 1 -type d)))  # find samples

for s in ${x[@]};  
do  
if [[ $(basename $s) != 'logs' ]]; 
then
    mkdir -p $data_path/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs
    echo "rsem-calculate-expression -p 8 --paired-end --bam --no-bam-output \
    $data_path/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/$(basename $s)/Aligned.toTranscriptome.out.bam \
    $pipelines_path/rsem/RSEM-1.3.3/files/$assembly/$annot/index \
    $data_path/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts" \
    | qsub -V -l nodes=1,mem=35gb,vmem=35gb,walltime=2:00:00 -j oe -d $data_path/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs/ \
    -N RSEM-bulk-T47D-$(basename $s);
fi




### + RSEM quantifications already launched by the bioinfo plateforme !!! to download on genomic plateform website (needs autentification) !!! :

###### DSP550 
###### DSP589 
###### DSP1111









########################################################################################
############################## INDEXING BAM FILES FOR IGV ##############################
########################################################################################


# sc T47D DSP762
################

samtools index $data_path/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt/Sample_T47D/Aligned.sortedByCoord.out.bam

# bulk T47D DSP280 (NI_E2)
##########################

samtools index $data_path/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/NI_E2/Aligned.sortedByCoord.out.bam


# bulk T47D DSP550 (T47D-pBABE-n1)
##################################

samtools index $data_path/iric/bulk/dsp550/downloaded-align/star/$assembly/$annot/raw-reads/T47D-pBABE-n1/Aligned.sortedByCoord.out.bam

# bulk ZR75 DSP550 (ZR75-pBABE-n1)
##################################

samtools index $data_path/iric/bulk/dsp550/downloaded-align/star/$assembly/$annot/raw-reads/ZR75-pBABE-n1/Aligned.sortedByCoord.out.bam


# bulk T47D DSP1111 (T-47D_pMIG_N1)
#############################

samtools index $data_path/iric/bulk/dsp1111/downloaded-align/star/$assembly/$annot/raw-reads/T-47D_pMIG_N1/Aligned.sortedByCoord.out.bam
