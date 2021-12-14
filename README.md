# RTCpredictor
RTCpredictor is a software tool for the prediction of read-through chimeric RNAs from RNA-Seq data.

**System Requirements:**
Linux or Mac OS with at least 5GB RAM. Also works on Windows OS via Cygwin or WSL (Windows Subsystem for Linux).

**Installation Instructions:**
This software requires 4 dependencies which needs to be installed/fulfilled prior to running it.

**Dependencies:**
<br>**1. ripgrep (version >=12.1.1)**
<br>Download pre-installed binary file of ripgrep software from https://github.com/BurntSushi/ripgrep as per your operating system (Linux/Windows/MacOS). Extract the compressed file and make sure the script 'rg' is in your PATH before running RTCpredictor software. Let's assume the full path of 'rg' binary script is /home/username/software/ripgrep/rg, then you can add it to the PATH by using following command:
<br>export PATH=/home/username/software/ripgrep:$PATH
<br>You need to use the above command every time you open a new terminal. In order to avoid it, add the above line to your ~/.bashrc file. An easy way to install ripgrep on MacOS is by using brew (command: brew install ripgrep) and on ubuntu >=18.10 (Linux) is by using apt-get (command: apt-get install ripgrep)

<br>**2. perl (>=v5.16.3).**
<br>perl is by default available in linux and macOS systems and needs to be installed in Cygwin in windows OS. RTCpredictor has not been tested using perl versions <v5.16.3 but it is expected to run on older versions of perl as well.

<br>**3. Parallel::ForkManager (version >=2.02)**
<br>Install this perl module. If you don't have admin privileges, we have provided an easy installation script (installParallelForkManager.bash) to install this perl module locally. Just run the script once using following command and you should be good to go:
<br>bash installParallelForkManager.bash

<br>**4. Database files.**
<br>We have provided database files for hg19 and hg38 genomes. [Click here](https://zenodo.org/record/5663811/files/hg19_data.tgz) to download hg19 database files (file size 116.2 MB). [Click here](https://zenodo.org/record/5663811/files/hg38_data.tgz) to download hg38 database files (file size 130.7 MB).

**Running the software:**
<br>Running the following command will provide all the options to run the software:
<br><br>perl RTCpredictor.pl -h
<br><br>Additionally, test fastq files (single and paired), small size database files from hg19 (for fast testing purpose) are also provided along with the expected result files. [Click here](https://zenodo.org/record/5663811/files/test_files.tgz) to download test files (file size 3.8 KB).
