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
<br>**(c) annotation file named 'genes_transcripts_exons.tsv'**: Detailed steps for making this file for hg38 genome as an example, is given below. The same steps can be followed for any other genome of interest.
<br>**(i)** Go to Ensembl website http://useast.ensembl.org/index.html
<br>**(ii)** Start a new BioMart instance by clicking on 'New' on the top left side of the webpage.
<br>**(iii)** Click on BioMart on the top menu
<br>**(iv)** Click on dropdown menu 'CHOOSE DATABASE'. Select 'Ensembl Genes 105'
<br>**(v)** Click on dropdown menu 'CHOOSE DATASET'. Select 'Human genes (GRCh38.p13)'
<br>**(vi)** On the left side menu, click on 'Attributes'. Among the options appearing on the top in bold, select 'Structures'.
<br>**(vii)** Now expand both sections 'GENE:' and 'EXON:' by clicking on '+' sign. Unselect all the checkboxes in both sections.
<br>**(viii)** Now select following checkboxes from section 'GENE:' in the same order as mentioned here: (1) 'Chromosome/scaffold name', (2) 'Strand', (3) 'Gene stable ID', (4) 'Gene name', (5) 'Gene start (bp)', (6) 'Gene end (bp)', (7) 'Transcript stable ID', (8) 'Transcript start (bp)', (9) 'Transcript end (bp)'.
<br>**(ix)** Now select following checkboxes from section 'EXON:' in the same order as mentioned here: (10) 'Exon stable ID', (11) 'Exon region start (bp)', (12) 'Exon region end (bp)'.
<br>**(x)** Please note that the order of selection of the checkboxes is very important. Please follow the order from number 1 to 12.
<br>**(xi)** On the top left side of the web page, click on 'Results'. In the resulting page, when asked for 'Export all results to', select 'File' and 'TSV' in the dropdown menu and select 'Unique results only' checkbox.
<br>**(xii)** Now click on 'Go'. A file named 'mart_export.txt' will be downloaded on your computer.
<br>**(xiii)** Rename this file to 'genes_transcripts_exons.tsv'
<br>**(d) paralog file named 'paralogs.tsv'**: Detailed steps for making this file for hg38 genome as an example, is given below. The same steps can be followed for any other genome of interest.
<br>**(i)** Go to Ensembl website http://useast.ensembl.org/index.html
<br>**(ii)** Click on BioMart on the top menu
<br>**(iii)** Click on dropdown menu 'CHOOSE DATABASE'. Select 'Ensembl Genes 105'
<br>**(iv)** Click on dropdown menu 'CHOOSE DATASET'. Select 'Human genes (GRCh38.p13)'
<br>**(v)** On the left side menu, click on 'Attributes'. Among the options appearing on the top in bold, select 'Homologues (Max select 6 orthologues)'.
<br>**(vi)** Now expand section 'GENE:' by clicking on '+' sign. Only 'Gene stable ID' checkbox should be selected and unselect all other checkboxes.
<br>**(vii)** Now expand section 'PARALOGUES:' by clicking on '+' sign and select 'Mouse paralogue gene stable ID' checkbox.
<br>**(viii)** On the top left side of the web page, click on 'Results'. In the resulting page, when asked for 'Export all results to', select 'File' and 'TSV' in the dropdown menu and select 'Unique results only' checkbox.
<br>**(ix)** Now click on 'Go'. A file named 'mart_export.txt' will be downloaded on your computer. Rename this file to 'paralogs.tsv'
<br>**(e) genome fasta file named 'genome.fa'**: Genome fasta file needs to be downloaded from ensembl. For hg38 genome, the file can be downloaded from http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz. After downloading, rename this file as 'genome.fa'. Moreover, the file needs to be processed in such a way that it's header should only have chromosome number and therefore the extra information in the header ( like 'dna:chromosome...') should be deleted. This can be done by running following command in the terminal
<br>perl -pi -e 's/dna:.+\n/\n/g' genome.fa
<br>**(f) (Optional) dbSNP bed file for the genome of interest**. For hg38 genome, the file can be downloaded from https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/. The individual files for each chromosome need to be concatenated/joined together to make a single bed file followed by sorting the file by coordinates. This can be done by running following commands in the terminal
<br>gunzip -c bed_chr_1.bed.gz bed_chr_2.bed.gz bed_chr_3.bed.gz bed_chr_4.bed.gz bed_chr_5.bed.gz bed_chr_6.bed.gz bed_chr_7.bed.gz bed_chr_8.bed.gz bed_chr_9.bed.gz bed_chr_10.bed.gz bed_chr_11.bed.gz bed_chr_12.bed.gz bed_chr_13.bed.gz bed_chr_14.bed.gz bed_chr_15.bed.gz bed_chr_16.bed.gz bed_chr_17.bed.gz bed_chr_18.bed.gz bed_chr_19.bed.gz bed_chr_20.bed.gz bed_chr_21.bed.gz bed_chr_22.bed.gz bed_chr_X.bed.gz bed_chr_Y.bed.gz > tmp_dbSNP.bed
<br>sort -k1,1 -k2,2n tmp_dbSNP.bed > dbSNP.bed
<br>**(g) (Optional) VCF file for the genome of interest**. For hg38 genome, the file can be downloaded from https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz. Please unzip the file using following command on the terminal
<br>gunzip 00-All.vcf.gz
**Running the script :**
<br>Running the following command will provide all the options to run the script:
<br><br>perl make_RTCpredictor_db.pl -h
