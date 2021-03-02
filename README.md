# Deciphering transcriptional networks underlying breast cancer subtype specification using single cell RNA-seq.

For all the sub-projects, you will need first to install all the pipelines and genome / transcriptome references required. To do so, you must run from this root repository the following :

./INSTALL.sh

Each subproject contains all the processing steps in a folder called 'analysis'. E.g., the first step 'count-matrix-generation' will download fastq files, reorder them, perform trimming of the reads and quantification, storing all the data in a folder called 'data'.

All folders containg heavy files, such as references, pipelines or data folders are not gonna be pushed with git (see .gitignore files).
