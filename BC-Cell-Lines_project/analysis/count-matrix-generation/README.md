## Generation of the count matrices

<br />

* To generate the count matrices (for both DSP762 and DSP992), run :

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; nohup ./fastq-to-mtx_DSP762.sh &> ./nohup_fastq-to-mtx_762.out < /dev/null &    <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; nohup ./fastq-to-mtx_DSP992.sh &> ./nohup_fastq-to-mtx_992.out < /dev/null &

&nbsp;&nbsp; (nohup allows you to prevent commands from being terminated when you log out or exit the terminal)

<br />

* If you have a compute cluster with torque installed, you can instead use qsub instead by running the following commands :

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ./fastq-to-mtx_DSP762.sh 1 <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ./fastq-to-mtx_DSP992.sh 1
