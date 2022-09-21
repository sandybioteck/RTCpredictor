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
<br>We have provided database files for hg19 and hg38 genomes. [Click here](https://zenodo.org/record/6407111/files/hg19_data.tgz) to download hg19 database files (file size 116.2 MB). [Click here](https://zenodo.org/record/6407111/files/hg38_data.tgz) to download hg38 database files (file size 130.7 MB).

**Running the software:**
<br>Running the following command will provide all the options to run the software:
<br><br>perl RTCpredictor.pl -h
<br><br>Additionally, test fastq files (single and paired), small size database files from hg19 (for fast testing purpose) are also provided along with the expected result files. [Click here](https://zenodo.org/record/6407111/files/test_files.tgz) to download test files (file size 3.8 KB).

**Making database files for any genome of interest**
<br>perl script 'make_RTCpredictor_db.pl' can be used to make the database files for any genome of interest. However, this script has some dependencies which needs to be installed/fulfilled prior to running it. Following are the dependencies of the script:
<br>**(a) blastn (version 2.7.1+ or later)**: pre-installed binary file of blastn compatible with any operating system (Linux/Windows/MacOS) can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
<br>**(b) bedtools (version 2.27.1 or later)**: pre-installed binary file of bedtools compatible with any operating system (Linux/Windows/MacOS) can be downloaded from https://github.com/arq5x/bedtools2
<br>**(c) annotation file named 'genes_transcripts_exons.tsv'**: Detailed steps for making this file is given below:
<br>Go to Ensembl website http://useast.ensembl.org/index.html
<br>Start a new BioMart instance by clicking on 'New' on the top left side of the webpage.
<br>Click on BioMart on the top menu
<br>Click on dropdown menu 'CHOOSE DATABASE'. Select 'Ensembl Genes 105'
<br>Click on dropdown menu 'CHOOSE DATASET'. Select 'Human genes (GRCh38.p13)'
<br>On the left side menu, click on 'Attributes'. Among the options appearing on the top in bold, select 'Structures'. Now expand both sections 'GENE:' and 'EXON:' by clicking on '+' sign. Unselect all the checkboxes in both sections. Now select following checkboxes from section 'GENE:' in the same order as mentioned here: (1) 'Chromosome/scaffold name', (2) 'Strand', (3) 'Gene stable ID', (4) 'Gene name', (5) 'Gene start (bp)', (6) 'Gene end (bp)', (7) 'Transcript stable ID', (8) 'Transcript start (bp)', (9) 'Transcript end (bp)'. Now select following checkboxes from section 'EXON:' in the same order as mentioned here: (10) 'Exon stable ID', (11) 'Exon region start (bp)', (12) 'Exon region end (bp)'. Please note that the order of selection of the checkboxes is very important. Please follow the order from number 1 to 12.
On the top left side of the web page, click on 'Results'. In the resulting page, when asked for 'Export all results to', select 'File' and 'TSV' in the dropdown menu and select 'Unique results only' checkbox. Now click on 'Go'. A file named 'mart_export.txt' will be downloaded on your computer. Rename this file to 'genes_transcripts_exons.tsv'
