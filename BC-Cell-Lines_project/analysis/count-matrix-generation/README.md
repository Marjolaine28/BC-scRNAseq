## Generation of the count matrices


To generate the count matrices (for both DSP762 and DSP992), run :


nohup ./fastq-to-mtx_DSP762.sh &> ./nohup_fastq-to-mtx_762.out < /dev/null &
nohup ./fastq-to-mtx_DSP992.sh &> ./nohup_fastq-to-mtx_992.out < /dev/null &

Nohup allows you to prevent commands from being terminated when you log out or exit the terminal.




If you have a compute cluster with torque installed, you can instead use qsub instead by running the following commands :

./fastq-to-mtx_DSP762.sh 1
./fastq-to-mtx_DSP992.sh 1
