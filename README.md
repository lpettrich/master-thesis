This is the Logfile for the Master thesis of Laura Pettrich 

Topic: Population History	

Duration: 03.01.2022 - 03.07.2022			       


# 1. First look at Crip4.0 Genome


seqkit stats [].fasta to get basic information about format type and length
  ### install seqkit
    i. Download latest tar.gz via firefox (link: https://github.com/shenwei356/seqkit/releases/tag/v2.1.0)
    
    ii. Decompress with tar -zxvf *.tar.gz

### samtools faidx  
   Example: 
    samtools faidx Chironomus_riparius_Mt.fasta Harlequin_PGI_Mt:10045-10833 > Chironomus_riparius_Mt.cox1.fasta 
    to extract COI region of the mitochondrial DNA of Crip  
    With the COI region your are able to identify the midge species more accurately. 
    Ann-Marie creates neighbour-joining trees and networks (using Geneious) with reference sequences.


## 1.1 Bloobtoolkit 
installation of Miniconda on cheops0 not possible, because only python 3.6 available on cheops and 3.8 at least needed
work on cheops1 for installing bloobtoolkit, module miniconda already available on cheops 1
NEEDS TO BE INSTALLED ON SCRATCH -> to much space needed!

    ssh -X lpettric@cheops1.rrz.uni-koeln.de

    module purge
    module load miniconda
    eval "$(conda shell.bash hook)"
--> Installation with conda did not work!!!

### Solution to Blobtools Viewer 


Don't install Blobtools locally!
Peter Heger installed it on Singularity as singularity image 
Path: /opt/rrzk/software/singularity_images/blobtoolkit_latest.sif

You need to mount your input and output directory with the parameter "-B"


#### Step 1 

i) Login on cheops1

ii) Modify the blobtools command so it fits to Singularity (i.e. "singularity exec") and mount your paths 

singularity exec -B /projects/ag-waldvogel/genomes/CRIP/crip4.0/:/schluppsi -B /projects/ag-waldvogel/genomes/CRIP/crip4.0/btk_output/:/outolino /opt/rrzk/software/singularity_images/blobtoolkit_latest.sif blobtools view --remote --out /outolino --view cumulative /schluppsi/blobtools


#### Step 2 

    i) Open new terminal

    ii) Login with information from "For remote access use ..." from the first step (username@remote_host = lpettric@cheops1.rrz.uni-koeln.de) 

    iii) Login with your password


#### Step 3 

    i) Copy link of "View dataset at .." from the first step

    ii) Paste it to your browser and load the webpage

 

## 1.2 Assembly stats

Use Cheops1
#### install assembly-stats
   /bin
   module load git
   git clone https://github.com/rjchallis/assembly-stats

#### Create a fresh conda env and install all needed packages (see: # See: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html and https://anaconda.org/

   module load miniconda
   eval "$(conda shell.bash hook)"
   conda create --prefix /scratch/lpettric/asst_env
   conda activate /scratch/lpettric/asst_env
   conda install python=3.8.5
   conda install perl
   conda install -c bioconda perl-json
   conda install -c bioconda perl-app-cpanminus
   conda install -c conda-forge -c bioconda busco=5.2.2

#### Create json file
   cd /projects/ag-waldvogel/genome/CRIP/crip4.0/
   perl /home/lpettric/bin/assembly-stats/pl/asm2stats.pl Chironomus_riparius_genome_010921.fasta > ./assembly-stats/Chironomus_riparius_genome_010921.assembly-stats.json
   perl /home/lpettric/bin/assembly-stats/pl/asm2stats.minmaxgc.pl Chironomus_riparius_genome_010921.fasta > ./assembly-stats/Chironomus_riparius_genome_010921.minmaxgc.assembly-stats.json

#### Add BUSCO (Rothamsted compared it to Insecta odb9)
C: Complete
S: Single-copy => not needed for assembly stats! "S":88.9
D: Duplicated
F: Fragmented
M: Missing
n: no. of genes => nicht sicher, ob eingesetze Zahl richtig
	"busco": { "C":96.9,
	            "D":8.0,
	            "F":0.1,
	            "M":3.0,
                    "n": 15212 },

 
#### Move to MinION computer:  scp lpettric@cheops.rrz.uni-koeln.de:/projects/ag-waldvogel/genomes/CRIP/crip4.0/assembly-stats/*.assembly-stats.json ./

json files created on Cheops using scripts as mentioned in instructions
installed assembly stats on MinION computer
integrated json files into html-script (also busco code)
open html-script with firefox on MinION computer --> It's okay if it's not a localhost website but opened as file!



## 1.3 R-package chromoMap

### install package on cheops0
   module purge
   module load R/4.0.2_gnu_mkl

### open R
   install.packages("chromoMap")





# 2. Pre-processing

Pre-processing of genome-scans of CRIP 

## Download reads from ENA (Project number: PRJEB24868)
Using FTP Client

### Use Cheops1

Command-line FTP clients allow you to interactively explore the FTP server and download data to your local computer. When asked for a username, use ‘anonymous’. When asked for a password, press the enter key to skip this.

   ftp ftp.sra.ebi.ac.uk
   Name: anonymous
   Password:
   ftp> cd vol1/run/ERR252/ERR2528543/
   ftp> get MF1_R1.paired.fastq_true.gz

In the above example, the ‘cd’ command is used to ‘change directory’ to the required directory. Then, the ‘get’ command is used to specify the file of interest. At any time, you can use ‘ls’ to view the content of the current directory. The command ‘pwd’ can be used to identify what the current directory is.


## 2.1 FastQC and MultiQC (0_fastqc.sh & 0_multiqc.sh)
### FastQC
#### Script example
   module purge
   module load openjdk/1.8.0_202

   cd /projects/ag-waldvogel/pophistory/CRIP/

   /home/lpettric/bin/FastQC/fastqc --threads 10 -o ./MF1/ ./MF1/MF1_R1.paired.fastq_true.gz ./MF1/MF1_R2.paired.fastq_true.gz



### MultiQC
#### Use Cheops1
#### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/multiqc_env/
   conda install multiqc

#### Script
   module purge

   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda activate /scratch/lpettric/conda/multiqc_env/

   cd /projects/ag-waldvogel/pophistory/CRIP/

   multiqc --filename multiqc_report_all_trimmed.html -o /projects/ag-waldvogel/pophistory/CRIP/ --file-list /projects/ag-waldvogel/pophistory/CRIP/multiqc-file.list

#### multiqc-file-list contains:
   ./MF1
   ./MF2
   ./MF3
   ./MF4
   ./MG2
   ./MG3
   ./MG4
   ./MG5
   ./NMF1
   ./NMF2
   ./NMF3
   ./NMF4
   ./SI1
   ./SI2
   ./SI3
   ./SI4
   ./SS1
   ./SS2
   ./SS3
   ./SS4


## 2.2 bwa mem (1_bwa-mem-script.sh)
### Script example
   module purge

   cd /projects/ag-waldvogel/pophistory/CRIP/

   /home/lpettric/bin/bwa/bwa mem -t 20 -M -R '@RG\tID:MF\tSM:MF1\tPL:ILLUMINA' /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta ./MF1/MF1_R1.paired.fastq_true.gz ./MF1/MF1_R2.paired.fastq_true.gz > ./bam-files/MF1_bwamem.sam

## 2.3 SAM to BAM (2_samtools-sort_modified.sh & 3_samtools-flagstat.sh)

### Script 1 - samtools sort
   module purge
   module load samtools

   cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

   ls -1 *_bwamem.sam | sed 's/_bwamem.sam//g' > list-crip
   while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-crip
   cat list-crip | /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools sort -@ 4 -o {}.bwamem.sort.bam {}.bam'

### Script 2 -samtools flagstat
   module purge
   module load samtools

   cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

   /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools flagstat {}.bwamem.sort.bam > {}.bwamem.sort.flagstat' < list-crip



## 2.4 Remove low quality alingment (4_samtools-remove-low-quality.sh)

### Script
   module purge
   module load samtools

   cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

   /home/lpettric/bin/parallel/bin/parallel -j 20 'samtools view -q 30 -f 0x0002 -F 0x0004 -F 0x0008 -b -o {}.bwamem.sort.q30.bam {}.bwamem.sort.bam' < list-crip



## 2.5 Remove duplicates (5_picard-remove-duplicates_newsyntax.sh

### Use Cheops1
### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/picard_env
   conda install -c bioconda picard 

### Script -> DON'T USE PARALLEL
   module purge
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda activate /scratch/lpettric/conda/picard_env

   cd /projects/ag-waldvogel/pophistory/CRIP/bam-files/

   picard MarkDuplicates -I MF1.bwamem.sort.q30.bam -O MF1.bwamem.sort.q30.rmd.bam -M MF1.bwamem.sort.q30.rmd.stat -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true


## 2.6 Collect mapping statistics

### Use Cheops1
### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/qualimap_env
   conda install -c bioconda qualimap 
   conda install R
   conda install samtools
   conda update qualimap

### Script
   module purge
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda activate /scratch/lpettric/conda/qualimap_env

   cd /projects/ag-waldvogel/pophistory/CRIP/bam-files/

   qualimap multi-bamqc -d ./list-qualimap -r --java-mem-size=4G -outformat PDF:HTML -outdir ./qualimap-output

### Content of list-qualimap
   MF1	./MF1.bwamem.sort.q30.rmd.bam	MF
   MF2	./MF2.bwamem.sort.q30.rmd.bam	MF
   MF3	./MF3.bwamem.sort.q30.rmd.bam	MF
   MF4	./MF4.bwamem.sort.q30.rmd.bam	MF
   ...



# 3. GATK

### Use Cheops1
### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/gatk_env
   conda install -c bioconda gatk4




# 4. Shapeit

### Use Cheops1
### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/shapeit_env
   conda install -c bioconda shapeit4=4.2.1

Core dumped!!!! calling shapeit4 does not work!!!
INSTEAD INSTALLED shapeit (v.2)
   conda install -c dranew shapeit 



# 5. Populations history anaylsis

## 5.1 iSMC (wrong package, instead eSMC)

### install iSMC
    i) Download latest tar.gz via firefox (link: https://github.com/gvbarroso/iSMC/releases/tag/v0.0.23)
    ii) Decompress with tar -zxvf *.tar.gz
    iii) Download additionals file using 
         git clone https://github.com/gvbarroso/iSMC


## 5.2 eSMC 

### Download github rep
   cd bin/
   git clone https://github.com/TPPSellinger/eSMC2.git
   
### Use Cheops1
### Create fresh conda env
   module load miniconda/py38_4.9.2
   eval "$(conda shell.bash hook)"
   conda create -p /scratch/lpettric/conda/eSMC_env
   conda install R
   conda install -c r r-devtools

### open R in eSMC conda env (eSMC_env)
  install.packages("rlang")
  install.packages("devtools")
  sessionInfo()


### open in R with eSMC conda env (eSMC_env)
	library(devtools)
	path="/home/lpettric/bin/eSMC2/eSMC2_1.1.5.tar.gz" # Path to the dowloaded eSMC package
	devtools::install_local(path)

   
