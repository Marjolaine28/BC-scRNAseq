## Installation

To install the pipelines versions used in this project, run the following command :

nohup ./INSTALL.sh ./logs/nohup_INSTALL.out < /dev/null &

Nohup prevents your commands from being terminated when you log out or exit the terminal ; the output of the command is redireted in a .out log file.


NOTE 1 : For any command you want to run, make sure you execute it from this github root repository.
NOTE 2 : All folders containg heavy files, such as references, pipelines or data folders are not gonna be pushed with git (see .gitignore files).



## Get data

To get the count matrices I used for the exploration and analysis of the scRNA-seq datasets (final workflow), you can run the following command :

nohup ./get-data_sc-workflow.sh __torque__ __key__ ./logs/nohup_get-data.out < /dev/null &

__Key__ is the identification that allows you acces to the data, you can find it in the URL that the genomic platform provides you to download your data (something like "https://genomique.iric.ca/FastQList?key=__the-key-you-want__&projectID=XXX&wget=1". __Torque__ is a ressource manager, if you use a compute cluster on witch it is installed, pass the value 1 otherwise pass the value 0. For example :

nohup ./get-data.sh __1__ __thisismykey__ ./logs/nohup_get-data.out < /dev/null &


To get the count matrices I used for Alevin parameters comparison (different k-mer sizes, etc.), you can run the following command :

nohup ./get-data_sc-workflow.sh __torque__ __key__ ./logs/nohup_get-data.out < /dev/null &


To get the count matrices (scRNA-seq as well as bulk RNA-seq datasets) I used for sanity check of the samples (IGV, correlations ...), you can run the following command :

nohup ./get-data_sc-workflow.sh __torque__ __key__ ./logs/nohup_get-data.out < /dev/null &


