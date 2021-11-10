# RTpred
RTpred is a software tool for the prediction of read-through chimeric RNAs from RNA-Seq data.

System Requirements:
Linux or Mac OS with at least 5GB RAM. Should work on Windows OS via Cygwin but didn't tested on it.

Installation Instructions:
This software requires 4 dependencies which needs to be installed prior to running it.
Dependencies:
1. ripgrep (version >=12.1.1)
Download pre-installed binary file of ripgrep software from https://github.com/BurntSushi/ripgrep. For MacOS, an easy way to install ripgrep is by using brew. (command: brew install ripgrep)
Extract the tar.gz file and make sure the script 'rg' is in your PATH before running RTpred software. Let's assume the full path of 'rg' binary script is /home/username/software/ripgrep/rg, then you can add it to the PATH by using following command:
export PATH=/home/username/software/ripgrep:$PATH
You need to use the above command every time you open a new terminal. In order to avoid it, add the above line to your ~/.bashrc file.

2. perl (>=v5.16.3). perl is by default available in linux and macOS systems. RTpred has not been tested using perl versions <v5.16.3 but it is expected to run on older versions of perl as well.

3. Parallel::ForkManager (version >=2.02)
Install this perl module. If you don't have admin privileges, we have provided an easy installation script (installParallelForkManager.bash) to install this perl module locally. Just run the script once using following command and you should be good to go:
bash installParallelForkManager.bash

4. database files. 

Running the software:
Running the following command will provide all the options to run the software:
perl RTpred.pl -h
Additionally, test fastq files (single and paired), small size database and exon files from hg19 (for fast testing purpose) are also provided along with the expected result files.
