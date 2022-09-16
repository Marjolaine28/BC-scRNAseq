# Workflow

## Installation

To install the pipelines versions used in this project, run the following command :

nohup ./INSTALL.sh ./logs/nohup_INSTALL.out < /dev/null &

Nohup prevents your commands from being terminated when you log out or exit the terminal ; the output of the command is redireted in a .out log file.


NOTE 1 : For any command you want to run, make sure you execute it from this github root repository.
NOTE 2 : All folders containg heavy files, such as references, pipelines or data folders are not gonna be pushed with git (see .gitignore files).



## Preproccessing

### Get count matrices

#### All steps in one script

To get the count matrices I used for the exploration and analysis presented in the notebooks, you can run the following command :

nohup ./get-data.sh __torque__ __key__ ./logs/nohup_get-data.out < /dev/null &

__Key__ is the identification that allows you acces to the data, you can find it in the URL that the genomic platform provides you to download your data (something like "https://genomique.iric.ca/FastQList?key=__the-key-you-want__&projectID=XXX&wget=1". __Torque__ is a ressource manager, if you use a compute cluster on witch it is installed, pass the value 1 otherwise pass the value 0. For example :

nohup ./get-data.sh __1__ __thisismykey__ ./logs/nohup_get-data.out < /dev/null &


#### Detailled steps

The script get-data.sh menionned above wraps differents steps of the preprocessing workflow, that you might want to run separatly. Those steps are detailled below.

##### Get references
##### Download the raw FASTQ files
##### Trimming
##### FASTQC
##### Whitelisting
##### Run Alevin

### Normalization


### How to choose parameters ?


## Exploraory analysis


