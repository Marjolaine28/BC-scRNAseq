#!/bin/bash

torque=${1:-"1"}
path=${2:-"$(git root)"}
key=${3:-"2212-e9fcdc80cfd658683c786d89297a28fd"}







#################################################################################
######################### PREPARE WORKING ENVIRONNEMENT #########################
#################################################################################


export PATH="$path/scripts/bash/utils/:$PATH" 
export PATH="$path/scripts/bash/quantification/:$PATH"
export PATH="$path/pipelines/sratoolkit/sratoolkit.3.0.0-centos_linux64/bin/:$PATH"
export PATH="$path/pipelines/salmon/salmon-1.4.0/bin/:$PATH"
export PATH="$path/pipelines/star/2.7.10a/bin/:$PATH"
export PATH="$path/pipelines/samtools/samtools-1.16/bin/:$PATH"
export PATH="$path/pipelines/rsem/RSEM-1.3.3/:$PATH"



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


if [[ ! -d $path/data/iric/bulk/dsp280/downloaded-fastqs ]]; then
    wget --no-check-certificate -O - \
    "https://genomique.iric.ca/FastQList?key=${key}&projectID=280&wget=1" \
    | wget --no-check-certificate -P $path/data/iric/bulk/dsp280/downloaded-fastqs -cri -
fi
merge-fastq-iric.sh $path/data/iric/bulk/dsp280/downloaded-fastqs $path/data/iric/bulk/dsp280/merged-fastqs



# 56 cell lines #
#################

mkdir $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded

x=$(cat $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SRR_Acc_List.txt) ; 
for d in $x; do 
    prefetch $d ; 
    fasterq-dump --split-files $d; 
    rm -r $d; 
    gzip ${d}_*.fastq; 
done

sed -i -e 's/T47D Kbluc/T47D_Kbluc/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/MDAMB175/MDAMB175VII/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/SUM1315/SUM1315MO2/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt

x=($(echo $(cut -f1 -d',' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))
y=($(echo $(cut -f8 -d',' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))

for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then mkdir -p ${y[$i]} ; mv ${x[$i]}*1.fastq.gz ${y[$i]}/${x[$i]}*R1.fastq.gz ; mv ${x[$i]}*2.fastq.gz ${y[$i]}/${x[$i]}*R2.fastq.gz ; fi ; done
for i in "${!x[@]}"; do if [[ $i != 0 ]] ; then echo ${y[$i]} ; echo ${x[$i]} ; fi ; done


merge-fastq-iric.sh $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/raw-fastqs








#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################


if [[ ! -f $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz ]]; 
then 
    mkdir -p $path/data/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $path/data/references/$assembly
else
    echo "Reference already $assembly available."
fi


annot=annotation__gencode/gencode_34

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz ]]; 
then 
    mkdir -p $path/data/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $path/data/references/$assembly/$annot
else
    echo "Reference already $annot available."
fi


annot=annotation__gencode/gencode_37

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v37.transcripts.fa.gz ]]; 
then 
    mkdir -p $path/data/references
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v37.transcripts.fa.gz -P $path/data/references/$assembly/$annot
else
    echo "Reference already $annot available."
fi





##################################################################################
######################## GET TRANSCRIPTS TO GENES MAPPING ########################
##################################################################################


annot=annotation__gencode/gencode_34

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv ]]; then
    zcat $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/|/\t/g' $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
    sed -i 's/>//g' $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv
fi


annot=annotation__gencode/gencode_37

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv ]]; then
    zcat $path/data/references/$assembly/$annot/gencode.v37.transcripts.fa.gz | grep '>' | cut -d"|" -f1,2 > $path/data/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
    sed -i 's/|/\t/g' $path/data/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
    sed -i 's/>//g' $path/data/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv
fi






#############################################################################
############################## SALMON INDEXING ##############################    
#############################################################################


annot=annotation__gencode/gencode_34
output_folder_index=$path/pipelines/salmon/salmon-1.4.0/files/$assembly/$annot 



# Gets decoys with gencode 34 #
###############################

mkdir -p $output_folder_index/decoys

if [[ ! -f $output_folder_index/decoys/decoys.txt && ! -f $output_folder_index/decoys/gentrome.fa ]] 
then
    salmon-get-decoys.sh -g $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz \
    -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys
fi 



# Indexing with decoys & gencode 34 #
#####################################

if [[ -d $output_folder_index/decoys/index_k31 ]]; then
    echo "decoys/index_k31 already generated.
    
    "
else
    if [[ $torque = 0 ]]; then
        export PATH="$path/pipelines/salmon/salmon-1.4.0/bin:$PATH"
        salmon index -k 31 -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -i $output_folder_index/decoys/index_k31 --gencode

    else
        qsub-salmon-indexing.sh -w 5:00:00 -m 150gb -t $output_folder_index/decoys/gentrome.fa.gz -d $output_folder_index/decoys/decoys.txt -o $output_folder_index/decoys -k 31 -p $path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
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
        salmon index -k 31 -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -i $output_folder_index/no-decoys/index_k31 --gencode
    else
        qsub-salmon-indexing.sh -w 1:00:00 -m 10gb -t $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz -o $output_folder_index/decoys -k 31 -p $path/pipelines/salmon/salmon-1.4.0 -args '--gencode'
        pids_index_no_decoys_k31=$(qstat -u $USER | tail -n 1 | awk '{print $1}' | cut -d"." -f1)
    fi
fi







###########################################################################
############################## STAR INDEXING ##############################    
###########################################################################


# Create STAR index with GRCh38-hg38 assembly and gencode 37 annotation #
#########################################################################

annot=annotation__gencode/gencode_37
output_folder_index=$path/pipelines/star/2.7.10a/files/$assembly/$annot


echo "STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $output_folder_index/index \
--genomeFastaFiles $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile $path/data/references/$assembly/$annot/gencode.v37.annotation.gtf" \
| qsub -V -l nodes=1,mem=200gb,vmem=200gb,walltime=48:00:00 -j oe -d $path/pipelines/star/2.7.10a/files/$assembly/$annot/index/logs/ -N "STAR-index-hg38-v37"









###########################################################################
############################## RSEM INDEXING ##############################    
###########################################################################


# Create RSEM reference with Gencode 37
#######################################

annot=annotation__gencode/gencode_37
output_folder_index=$path/pipelines/rsem/RSEM-1.3.3/files/$assembly/$annot 

rsem-prepare-reference --gtf $path/data/references/$assembly/$annot/gencode.v37.annotation.gtf \
$path/data/references/$assembly/GRCh38.primary_assembly.genome.fa $output_folder_index/index/ ####### QSUB ? WAIT PID ??









#############################################################################################################
##############################  MAPPING BULK DATA WITH SALMON FOR CORRELATIONS ##############################
#############################################################################################################


# 56 cells lines - decoys k31 & gencode 34 #
############################################

annot=annotation__gencode/gencode_34
output_folder_index=$path/pipelines/salmon/salmon-1.4.0/files/$assembly/$annot/index_k31 

qsub-all-fastqs.sh -w 20:00:00 -m 100gb -h 0 -r $path/scripts/bash/quantification/run-salmon.sh -f $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/raw-fastqs \
-o $path/data/public/RNA-seq/bulk/BC-Cell-Lines-Panel_GSE48213/quant/salmon/$assembly/$annot/raw-reads/pc-decoys-k31-1.4.0/raw-counts -l "PE" -s "all" -p $path/pipelines/RNA-seq/salmon/salmon-1.4.0 \
-args "-l A -i $output_folder_index -g $path/data/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"


# MCF7 (DSP356) - no decoys k31 & gencode 34 #
##############################################


qsub-all-fastqs.sh -w 10:00:00 -m 60gb -r $path/scripts/bash/quantification/run-salmon.sh -f $path/data/iric/bulk/dsp356/raw-fastqs \
-o $path/data/iric/bulk/dsp356/quant/salmon/$assembly/$annot/raw-reads/no-decoys-k31-1.4.0/raw-counts -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
-p $path/pipelines/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/no-decoys/index_k31 \
-g $output_folder_index -g $path/data/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"


# MCF7 (DSP356) - decoys k31 & gencode 34 #
###########################################


qsub-all-fastqs.sh -w 10:00:00 -m 60gb -r $path/scripts/bash/quantification/run-salmon.sh -f $path/data/iric/bulk/dsp356/raw-fastqs \
-o $path/data/iric/bulk/dsp356/quant/salmon/$assembly/$annot/raw-reads/decoys-k31-1.4.0/raw-counts -l "PE" -s "Sample_E-1  Sample_E-2  Sample_E-3 Sample_V-1  Sample_V-2  Sample_V-3" \
-p $path/pipelines/salmon/salmon-1.4.0 -args "-l A -i $output_folder_index/decoys/index_k31 \
-g $output_folder_index -g $path/data/references/$assembly/$annot/gencode.v34.pc_transcripts_txp2gene.tsv --validateMappings"








####################################################################################################
############################## STAR ALIGNMENTS FOR CORRELATIONS & IGV ##############################
####################################################################################################


annot=annotation__gencode/gencode_37


# Align sc T47D (DSP762) in bulk mode with STAR (just for IGV, so no need to quantify after)
############################################################################################


qsub-all-fastqs.sh -w 30:00:00 -m 100gb -r $path/scripts/bash/quantification/run-star.sh -f $path/data/iric/sc/dsp762/trimmed-fastqs \ #### CHECK TRIMMED PATH
-o $path/data/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt -l "PEr2" -s "Sample_T47D" \
-args "--genomeDir $path/pipelines/star/2.7.10a/files/$assembly/$annot/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"





# Align bulk T47D (DSP280) with STAR
####################################


qsub-all-fastqs.sh -w 30:00:00 -m 100gb -h 0 -r $path/scripts/bash/quantification/run-star.sh -f $path/data/iric/bulk/dsp280/raw-fastqs \
-o $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads -l "PE" -s "all" \
-args "--genomeDir $path/pipelines/star/2.7.10a/files/$assembly/$annot/index/ \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard"



# Load from server in IGV : https://genomique.iric.ca/IGVRegistry?key=2212$e9fcdc80cfd658683c786d89297a28fd




####################################################################################################
############################## RSEM QUANTIFICATION FOR CORRELATIONS & IGV ##########################
####################################################################################################


annot=annotation__gencode/gencode_37


# RSEM quant for bulk T47D (DSP280)
###################################

x=($(echo $(find $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads -mindepth 1 -maxdepth 1 -type d)))  # find samples

for s in ${x[@]};  
do  
if [[ $(basename $s) != 'logs' ]]; 
then
    mkdir -p $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs
    echo "rsem-calculate-expression -p 8 --paired-end --bam --no-bam-output \
    $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/$(basename $s)/Aligned.toTranscriptome.out.bam \
    $path/pipelines/rsem/RSEM-1.3.3/files/$assembly/$annot/index \
    $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts" \
    | qsub -V -l nodes=1,mem=35gb,vmem=35gb,walltime=2:00:00 -j oe -d $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs/ \
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

samtools index $path/data/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt/Sample_T47D/Aligned.sortedByCoord.out.bam

# bulk T47D DSP280 (NI_E2)
##########################

samtools index $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/NI_E2/Aligned.sortedByCoord.out.bam


# bulk T47D DSP550 (T47D-pBABE-n1)
##################################

samtools index $path/data/iric/bulk/dsp550/downloaded-align/star/$assembly/$annot/raw-reads/T47D-pBABE-n1/Aligned.sortedByCoord.out.bam

# bulk ZR75 DSP550 (ZR75-pBABE-n1)
##################################

samtools index $path/data/iric/bulk/dsp550/downloaded-align/star/$assembly/$annot/raw-reads/ZR75-pBABE-n1/Aligned.sortedByCoord.out.bam


# bulk T47D DSP1111 (T-47D_pMIG_N1)
#############################

samtools index $path/data/iric/bulk/dsp1111/downloaded-align/star/$assembly/$annot/raw-reads/T-47D_pMIG_N1/Aligned.sortedByCoord.out.bam
