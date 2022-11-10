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

export PATH="$path/pipelines/samtools/samtools-1.16/bin/:$PATH"




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
# ################

mkdir -p $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs

x=$(cat $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/test_SRR_Acc_List.txt) ; 
for d in $x; do 
    prefetch $d -O $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs ; 
    fasterq-dump --split-files $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/$d \
        --outdir $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/$d -t $path/tmp ; 
    gzip $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/$d/${d}_*.fastq ; 
done

sed -i -e 's/T47D Kbluc/T47D_Kbluc/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/MDAMB175/MDAMB175VII/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt
sed -i -e 's/SUM1315/SUM1315MO2/g' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt

x=($(echo $(cut -f1 -d',' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))
y=($(echo $(cut -f8 -d',' $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/sra-files/SraRunTable.txt)))



for i in "${!x[@]}"; do 
    if [[ $i != 0 ]] ; then
        mv $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/${y[$i]}/${y[$i]}*1.fastq.gz \
            $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/${y[$i]}/${x[$i]}_R1.fastq.gz

        mv $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/${y[$i]}/${y[$i]}*2.fastq.gz \
            $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs/${y[$i]}/${x[$i]}_R2.fastq.gz
    fi
done


merge-fastq-iric.sh $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/downloaded-fastqs $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/merged-fastqs








#########################################################################################################
############################## DOWNLOAD REFERENCE GENOME AND TRANSCRIPTOME ##############################
#########################################################################################################


if [[ ! -f $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz ]]; 
then 
    mkdir -p $path/data/references/$assembly
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz -P $path/data/references/$assembly
else
    echo "Reference $assembly already available."
fi


annot=annotation__gencode/gencode_34

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v34.transcripts.fa.gz || ! -f $path/data/references/$assembly/$annot/gencode.v34.basic.annotation.gtf.gz ]]; 
then 
    mkdir -p $path/data/references/$assembly/$annot
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.basic.annotation.gtf.gz -P $path/data/references/$assembly/$annot
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz -P $path/data/references/$assembly/$annot
else
    echo "Reference $annot already available."
fi



annot=annotation__gencode/gencode_37

if [[ ! -f $path/data/references/$assembly/$annot/gencode.v37.transcripts.fa.gz || ! -f $path/data/references/$assembly/$annot/gencode.v37.basic.annotation.gtf.gz ]]; 
then 
    mkdir -p $path/data/references/$assembly/$annot
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.basic.annotation.gtf.gz -P $path/data/references/$assembly/$annot
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz -P $path/data/references/$assembly/$annot
else
    echo "Reference $annot already available."
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


if [[ ! -f $path/data/references/$assembly/$annot/gencode.v37.transcripts_gene2txp.tsv ]]; then
    awk '{print $2,$1}' $path/data/references/$assembly/$annot/gencode.v37.transcripts_txp2gene.tsv > $path/data/references/$assembly/$annot/gencode.v37.transcripts_gene2txp.tsv
fi



############################################################################
############################# SALMON INDEXING ##############################    
############################################################################


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









###########################################################################
############################## STAR INDEXING ##############################    
###########################################################################


# Create STAR index with GRCh38-hg38 assembly and gencode 37 annotation
#######################################################################

annot=annotation__gencode/gencode_37
output_folder_index=$path/pipelines/star/STAR-2.7.10a/files/$assembly/$annot

echo "Extracting references ..."
# gunzip $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa.gz;
# gunzip $path/data/references/$assembly/$annot/gencode.v37.basic.annotation.gtf.gz;


command="export PATH="$path/pipelines/star/STAR-2.7.10a/source/:$PATH";
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $output_folder_index/index \
--genomeFastaFiles $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile $path/data/references/$assembly/$annot/gencode.v37.basic.annotation.gtf"

if [[ -d $output_folder_index/index ]]; then
    echo "STAR index already generated.
    
    "
else
    if [[ $torque = 0 ]]; then
        $command
    else
        mkdir -p $output_folder_index/index/logs/
        pid_index_star=$(echo $command | qsub -V -l nodes=1,mem=150gb,vmem=150gb,walltime=48:00:00 -j oe -d $output_folder_index/index/logs/ -N "STAR-index-hg38-v37")
    fi
fi






###########################################################################
############################## RSEM INDEXING ##############################    
###########################################################################


# Create RSEM reference with Gencode 37
#######################################

annot=annotation__gencode/gencode_37
output_folder_index=$path/pipelines/rsem/RSEM-1.3.3/files/$assembly/$annot


command="export PATH="$path/pipelines/rsem/RSEM-1.3.3/:$PATH";
rsem-prepare-reference --gtf $path/data/references/$assembly/$annot/gencode.v37.basic.annotation.gtf -p 8 --star \
--star-path $path/pipelines/star/STAR-2.7.10a/source $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa $output_folder_index/index/"


if [[ -d $output_folder_index/index ]]; then
    echo "RSEM index already generated.
    
    "
else
    if [[ $torque = 0 ]]; then
        $command
    else
        mkdir -p $output_folder_index/index/logs
        pid_index_rsem=$(echo $command | qsub -V -l nodes=1,mem=200gb,vmem=200gb,walltime=48:00:00 -j oe -d $output_folder_index/index/logs/ -N "RSEM-index")
    fi
fi


echo "Compressing  references ..."

# gzip $path/data/references/$assembly/GRCh38.primary_assembly.genome.fa
# gzip $path/data/references/$assembly/$annot/gencode.v37.basic.annotation.gtf







#############################################################################################################
##############################  MAPPING BULK DATA WITH SALMON FOR CORRELATIONS ##############################
#############################################################################################################


# 56 cells lines - decoys k31 & gencode 34 #
############################################

annot=annotation__gencode/gencode_34
output_folder_index=$path/pipelines/salmon/salmon-1.4.0/files/$assembly/$annot/index_k31 

wait_pid=$pids_index_decoys_k31

$submit -w 20:00:00 -m 100gb -h "$wait_pid" -r $path/scripts/bash/quantification/run-salmon.sh -f $path/data/public/bulk/BC-Cell-Lines-Panel_GSE48213/merged-fastqs \
-o $path/data/public/RNA-seq/bulk/BC-Cell-Lines-Panel_GSE48213/quant/salmon/$assembly/$annot/raw-reads/decoys-k31-1.4.0/raw-counts -l "PE" -s "all" -p $path/pipelines/RNA-seq/salmon/salmon-1.4.0 \
-args "-l A -i $output_folder_index -g $path/data/references/$assembly/$annot/gencode.v34.transcripts_txp2gene.tsv --validateMappings"










####################################################################################################
############################## STAR ALIGNMENTS FOR CORRELATIONS & IGV ##############################
####################################################################################################


annot=annotation__gencode/gencode_37


# Align sc T47D (DSP762) in bulk mode with STAR (just for IGV, so no need to quantify after)
############################################################################################


annot=annotation__gencode/gencode_37
output_folder_index=$path/pipelines/star/STAR-2.7.10a/files/$assembly/$annot

wait_pid="$pid_index_star"
mkdir -p $path/data/iric/sc/dsp762/align/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/trimmed-reads-cutadapt/T47D/logs/

command="export PATH="$path/pipelines/star/STAR-2.7.10a/source:$PATH"; \
    STAR --genomeDir $output_folder_index/index/ --runThreadN 8 \
    --readFilesIn $path/data/iric/sc/dsp762/trimmed-fastqs/cutadapt/Sample_T47D/trimmed-Sample_T47D_R2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix $path/data/iric/sc/dsp762/align/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/trimmed-reads-cutadapt/T47D/T47D_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard"

if [[ -s $path/data/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt/T47D/T47D_Aligned.toTranscriptome.out.bam ]]; then
    echo "T47D_Aligned.toTranscriptome.out.bam already exists.
    
    "
elif [[ $torque = 0 ]]; then
    $command
elif [[ $wait_pid = "" ]]; then
    echo $command | qsub -V -l nodes=1,mem=60gb,vmem=60gb,walltime=24:00:00 -j oe \
    -d $path/data/iric/sc/dsp762/align/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/trimmed-reads-cutadapt/T47D/logs/ \
    -N "STAR-align-T47D"
else
    echo $command | qsub -V -l nodes=1,mem=60gb,vmem=60gb,walltime=24:00:00 -j oe \
    -d $path/data/iric/sc/dsp762/align/star/human/assembly__GRCh38-hg38/annotation__gencode/gencode_37/trimmed-reads-cutadapt/T47D/logs/ \
    -W depend=afterok:$wait_pid -N "STAR-align-T47D"
fi






####################################################################################################
############################## RSEM QUANTIFICATION FOR CORRELATIONS & IGV ##########################
####################################################################################################


annot=annotation__gencode/gencode_37


# RSEM quant for bulk T47D (DSP280)
###################################


wait_pid="$pid_index_rsem"


$submit -w 5:00:00 -m 60gb -h "$wait_pid" -r $path/scripts/bash/quantification/run-rsem.sh -f $path/data/iric/bulk/dsp280/merged-fastqs \
    -o $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/raw-counts \
    -l "PE" -s "all" -p $path/pipelines/rsem/RSEM-1.3.3/ -savepids $path/tmp/qsub_pids \
    -args "$path/pipelines/rsem/RSEM-1.3.3/files/$assembly/$annot/index/ --star --star-path $path/pipelines/star/STAR-2.7.10a/source --star-gzipped-read-file --output-genome-bam -p 8" # index path must be in first position in args




# #### A enlever

# # x=($(echo $(find $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads -mindepth 1 -maxdepth 1 -type d)))  # find samples

# # for s in ${x[@]};  
# # do  
# # if [[ $(basename $s) != 'logs' ]]; 
# # then
# #     mkdir -p $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs
# #     echo "rsem-calculate-expression -p 8 --paired-end --bam --no-bam-output \
# #     $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/$(basename $s)/Aligned.toTranscriptome.out.bam \
# #     $path/pipelines/rsem/RSEM-1.3.3/files/$assembly/$annot/index \
# #     $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts" \
# #     | qsub -V -l nodes=1,mem=35gb,vmem=35gb,walltime=2:00:00 -j oe -d $path/data/iric/bulk/dsp280/quant/rsem-star/$assembly/$annot/raw-reads/$(basename $s)/raw-counts/logs/ \
# #     -N RSEM-bulk-T47D-$(basename $s);
# # fi




### + RSEM quantifications already launched by the bioinfo plateforme !!! to download on genomic plateform website (needs autentification) !!! :

###### DSP550 
###### DSP589 
###### DSP1111









########################################################################################
############################## INDEXING BAM FILES FOR IGV ##############################
########################################################################################


# sc T47D DSP762
################

samtools index $path/data/iric/sc/dsp762/align/star/$assembly/$annot/trimmed-reads-cutadapt/T47D/T47D_Aligned.sortedByCoord.out.bam

# # bulk T47D DSP280 (NI_E2)
# ##########################

# samtools index $path/data/iric/bulk/dsp280/align/star/$assembly/$annot/raw-reads/NI_E2/Aligned.sortedByCoord.out.bam



# THEN
# - Add index files into IGV for DSP280 and DSP762
# - DSP550 and DSP1111 (proceessed by the genomic platform) can be loaded from server in IGV : https://genomique.iric.ca/IGVRegistry?key=2212$e9fcdc80cfd658683c786d89297a28fd
