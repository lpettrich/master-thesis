This is the Logfile for the Master's Thesis by Laura Pettrich 

Topic: Population History	

Duration: 03.01.2022 - 03.07.2022			       

# ***Chironomus riparius***
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

Use Cheops1

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

## 5.2 MSMC2 

### Download msmc-tools and msmc2 (git clone)

#### Download and install dmd
     curl https://dlang.org/install.sh | bash -s
     ~/dlang/install.sh update

Run `source ~/dlang/dmd-2.099.0/activate` in your shell to use dmd-2.099.0.

This will setup PATH, LIBRARY_PATH, LD_LIBRARY_PATH, DMD, DC, and PS1.

Run `deactivate` later on to restore your environment.

#### Download latest GSL from GNU mirror (https://www.gnu.org/software/gsl/#downloading)
How to install: https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/
     mkdir /home/lpettric/bin/gsl
     ./configure --prefix=/home/lpettric/bin/gsl
     make
     make check 
     make install

You can then install the program by typing make in the directory. The resulting executable will be in the build/release subdirectory.

You have to manually adjust the Makefile if you have the GSL in a non-trivial location on your platform. = modify to /home/lpettric/bin/gsl/lib/


**For MSMC2 you need phased haplotype data of single chromosomes as input**


### 5.2.1 Preparation


#### Get names of chromosomes


     grep '^>' Chironomus_riparius_genome_010921.fasta | sed 's/>//' > chromosome-nr.txt

I will only use the chromosomes not the scaffolds of the crip4.0 chromosme


#### Index bam-files


     while read f; do samtools index $f".bwamem.sort.q30.rmd.bam" > $f".bwamem.sort.q30.rmd.bam.bai" ; done < list-crip

MF2.bwamem.sort.q30.rmd.bam seems to be broken needs to be redone! -> Check!



### 5.2.2 Extract chromosome information + index

a) Extract chromosome

     samtools faidx Chironomus_riparius_genome_010921.fasta Chr1 > Chr1_Chironomus_riparius_genome_010921.fasta
     samtools faidx Chironomus_riparius_genome_010921.fasta Chr2 > Chr2_Chironomus_riparius_genome_010921.fasta
     samtools faidx Chironomus_riparius_genome_010921.fasta Chr3 > Chr3_Chironomus_riparius_genome_010921.fasta
     samtools faidx Chironomus_riparius_genome_010921.fasta Chr4 > Chr4_Chironomus_riparius_genome_010921.fasta


b) Index 

     /home/lpettric/bin/bwa/bwa index Chr1_Chironomus_riparius_genome_010921.fasta
     /home/lpettric/bin/bwa/bwa index Chr2_Chironomus_riparius_genome_010921.fasta
     /home/lpettric/bin/bwa/bwa index Chr3_Chironomus_riparius_genome_010921.fasta
     /home/lpettric/bin/bwa/bwa index Chr4_Chironomus_riparius_genome_010921.fasta


### 5.2.3 Create mappability mask per chromosome using SNPable

a) Extract overlapping 145mers subsequences as artificial reads from Chr1

I have the same data as Ann-Marie so I choose the same length

     /home/lpettric/bin/seqbility-20091110/splitfa Chr1_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
     /home/lpettric/bin/seqbility-20091110/splitfa Chr2_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
     /home/lpettric/bin/seqbility-20091110/splitfa Chr3_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000
     /home/lpettric/bin/seqbility-20091110/splitfa Chr4_Chironomus_riparius_genome_010921.fasta 145 | split -l 20000000

Move into new subdirectory

b) Map artificial reads back to chromosomes

     /home/lpettric/bin/bwa/bwa aln -R 1000000 -O 3 -E 3 ../Chr1_Chironomus_riparius_genome_010921.fasta xaa > xaa_bwaaln.sai
    ... repeat for all xa? of each chromosome

c) Convert sai to sam

     /home/lpettric/bin/bwa/bwa samse ../Chr1_Chironomus_riparius_genome_010921.fasta xaa_Chr1_bwaaln.sai xaa > xaa_Chr1_bwaaln.sam
    ... repeat for all xa? of each chromosome

    bgzip and pileup all mappings of one chromosme: Chr1-all_bwaaln.sam.gz

==> NOW I HAVE SEVERAL MAPPINGS FOR ONE CHROMOSME --> NOT IDEAL --> merge xa? files to one file and repeat steps b) + c)
 
     cat xaa xab xac xad xae xaf xag > Chr1_145splits.fa
     cat xaa xab xac xad xae xaf > Chr2_145splits.fa
     cat xaa xab xac xad xae xaf > Chr3_145splits.fa
     cat xaa xab > Chr4_145splits.fa

b)   /home/lpettric/bin/bwa/bwa aln -R 1000000 -O 3 -E 3 Chr1/Chr1_Chironomus_riparius_genome_010921.fasta Chr1/145mer-subsequences/Chr1_145splits.fa > Chr1/145mer-subsequences/Chr1_145splits_bwaaln.sai 

c)   /home/lpettric/bin/bwa/bwa samse Chr1/Chr1_Chironomus_riparius_genome_010921.fasta Chr1/145mer-subsequences/Chr1_145splits_bwaaln.sai Chr1/145mer-subsequences/Chr1_145splits.fa > Chr1_145splits_bwaaln.sam

    gzip Chr1_145splits_bwaaln.sam

d) Generate rawMask
     gzip -dc Chr1_145splits_bwaaln.sam.gz | /home/lpettric/bin/seqbility-20091110/gen_raw_mask.pl > rawMask_Chr1_145.fa

e) Generate the final mask
     gen_mask -l 145 -r 0.5 rawMask_Chr1_145.fa > mask_Chr1_145_50.fa

length: 145bp and stringency: 0.5


final-masks can be found in a seperate directory

f) convert final-masks to .bed using makeMappabilitMask.py

change paths of input and output in script


### 5.2.4 Phase reads

OPTIONAL: Indel realingment (GATK3 not on GATK4 (with freebayes not neccessary))

You must split the dataset by chromosomes prior to phasing

#### a) Variant calling (you could also use freebayes and filter with vcftool but I will use script from msmc-tools (uses samtools pileup and bcftools call)

a.1) Try script from msmc-tools (bamCaller.py)

Get coverage statistics per chromosome

     while read f; do samtools depth -r Chr1 $f".bwamem.sort.q30.rmd.bam" | awk '{sum += $3} END {print sum / NR}' > $f".Chr1.cov" ; done < list-crip
     # repeat for every chromosome

Run bamCaller.py

     samtools mpileup -q 30 -Q 20 -C 50 -u -r <chr> -f <ref.fa> <bam> | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py <mean_cov> <out_mask.bed.gz> | gzip -c > <out.vcf.gz>
    
samtools:

    # q = Minimum mapping quality for an alignment to be used
    
    # Q = Minimum base quality for a base to be considered.
    
    # C = Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. 
    
    # u = uncompressed
    
    # r = Only generate pileup in region. Requires the BAM files to be indexed. If used in conjunction with -l then considers the intersection of the two requests.
    
    # f = fasta-ref
    
 bcftools:
 
    # c = consensus-caller
 
    # V = skip-variants snps|indels

Get summary list of coverage per chromosome per bam

    awk '{print $0 "\t" FILENAME}' *.Chr1.cov > Chr1.summary


Final command: Instead of samtools use bcftools because mpielup migrated to it

    module purge
    module load samtools/1.13
    module load python/3.4.3

    cd /projects/ag-waldvogel/pophistory/CRIP/bam-files

    # Chr1
    while read -r x y; do bcftools mpileup -q 30 -Q 20 -C 50 -r Chr1 -f /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta $y".bwamem.sort.q30.rmd.bam" | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py $x $y"_mask.bed.gz" | gzip -c > $y".vcf.gz" ;done < ./mean-coverage-chromosome/Chr1.summary

--> changed script, so that it is directly bgzipped

    # Chr2
    while read -r x y; do bcftools mpileup -q 30 -Q 20 -C 50 -r Chr2 -f /projects/ag-waldvogel/genomes/CRIP/crip4.0/Chironomus_riparius_genome_010921.fasta $y".bwamem.sort.q30.rmd.bam" | bcftools call -c -V indels | /home/lpettric/bin/msmc-tools/bamCaller.py $x $y"_Chr2_mask.bed.gz" | bgzip -c > $y"_Chr2.vcf.gz" ;done < ./mean-coverage-chromosome/Chr2.summary


merge vcf.files and zip and index them (because no reference panel)

     bgzip *.vcf --> bcftools only index if they are bgzipped
     
normal gzip worng! need to change to bgzip
command for already existing files of Chr1: 

       for f in *vcf*; do zcat $f | bgzip -c > $f".bgz" ; done
       
check type of vcf with:

      htsfile file.vcf.gz
      
index

          for f in *vcf.gz; do bcftools index $f; done

you get csi index
 
      bcftools merge --print-header *.vcf.gz # to get header info to see if there are recurring headers
      bcftools merge -O z -o merged_Chr1.vcf.gz  *.vcf.gz       

#### b) No reference panel, that's why we merged the vcf, now we filter to only have monoallelic and biallelic SNPs  

     bcftools view -M 2 -O z -o merged_biallelic_Chr1.vcf.gz merged_Chr1.vcf.gz 
     bcftools index merged_biallelic_Chr1.vcf.gz

#### c) Phasing with SHAPEIT --> find recombination rate (https://academic.oup.com/g3journal/article/10/4/1151/6026161#235660423: ranges between 0.04 and 0.07 on Chr3 less, use 0.04 since it is on Chr3 less --> no parameter with recombination rate in shapeit4)

Shapeit4.2 was modified by Peter Heger to remove AVX2 dependency

start shapeit main run (need to load boost and samtools)

     /home/lpettric/bin/shapeit4/bin/shapeit4.2 -I merged_biallelic_Chr1.vcf.gz -O ./phased/merged_biallelic_Chr1_phased.vcf.gz --sequencing --region Chr1 --log  /phased/shapeit_Chr1.log

after phasing all individuals together, seperate them again

      while read f; do bcftools view -s $f -O z -o $f"_Chr1_phased.vcf.gz" merged_biallelic_Chr1_phased.vcf.gz ; done < ../../../bam-files/list-crip


vcf.gz files indexed

      for f in *.vcf.gz; do bcftools index $f; done

#### e) Correct for missed genotypes
These files now contain the pashed alleles for each individual. However, the pre-phased files might contain more information that should not be lost.

Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones

Use --force-samples because phased and unphase have same headers

      while read f; do bcftools merge --force-samples ../$f"_Chr1.vcf.gz" $f"_Chr1_phased.vcf.gz" | awk '
    	  BEGIN {OFS="\t"}
 	     $0 ~ /^##/ {print}
	     $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
	     $0 !~ /^#/ {
 	     if(substr($11, 1, 3) != "./.")
  	     $10 = $11
 	     print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
      }' | bcftools view -O z > $f"_Chr1_phased_merged.vcf.gz" ; done < ../../../bam-files/list-crip

### 5.2.5 Create input-files (multihetsep)

	while read a b c d; do /home/lpettric/bin/msmc-tools//generate_multihetsep.py --mask=../$a"_Chr1_mask.bed.gz" \
                          --mask=..$b"_Chr1_mask.bed.gz" \
                          --mask=..$c"_Chr1_mask.bed.gz" \
                          --mask=..$d"_Chr1_mask.bed.gz" \
                          --mask=/projects/ag-waldvogel/pophistory/CRIP/masking/final-mask/mask_Chr1_145_50.bed.gz \
                          $a"_Chr1_phased_merged.vcf.gz" $b"_Chr1_phased_merged.vcf.gz" \
                          $c"_Chr1_phased_merged.vcf.gz" $d"_Chr1_phased_merged.vcf.gz" > "multihetsep_"$a"_"$b"_"$c"_"$d"_Chr1.txt"; done < ../../../msmc2/list-populations

#### 5.2.6 First runs of msmc2 = one population only
USE MUTATION RATE FROM ANN-MARIE'S PAPER '21                  
4.27x10^-9 #from Waldvogel & Pfenninger 2021; mutation rate   
not for msmc run but to get real time data                    

##### Wiederholung mit anderem -p (run6)
On Cheops1
	module purge

	cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/

	## MG
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o MF_Chr1-4_msmc_output ./multihetsep-Chr1/multihetsep_MF1_MF2_MF3_MF4_Chr1.txt ./multihetsep-Chr2/multihetsep_MF1_MF2_MF3_MF4_Chr2.txt ./multihetsep-Chr3/multihetsep_MF1_MF2_MF3_MF4_Chr3.txt ./multihetsep-Chr4/multihetsep_MF1_MF2_MF3_MF4_Chr4.txt &
	wait
	## MG
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o MG_Chr1-4_msmc_output ./multihetsep-Chr1/multihetsep_MG2_MG3_MG4_MG5_Chr1.txt ./multihetsep-Chr2/multihetsep_MG2_MG3_MG4_MG5_Chr2.txt ./multihetsep-Chr3/multihetsep_MG2_MG3_MG4_MG5_Chr3.txt ./multihetsep-Chr4/multihetsep_MG2_MG3_MG4_MG5_Chr4.txt &
	wait
	## NMF
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o NMF_Chr1-4_msmc_output ./multihetsep-Chr1/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr1.txt ./multihetsep-Chr2/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr2.txt ./multihetsep-Chr3/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr3.txt ./multihetsep-Chr4/multihetsep_NMF1_NMF2_NMF3_NMF4_Chr4.txt &
	wait
	## SI
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o SI_Chr1-4_msmc_output ./multihetsep-Chr1/multihetsep_SI1_SI2_SI3_SI4_Chr1.txt ./multihetsep-Chr2/multihetsep_SI1_SI2_SI3_SI4_Chr2.txt ./multihetsep-Chr3/multihetsep_SI1_SI2_SI3_SI4_Chr3.txt ./multihetsep-Chr4/multihetsep_SI1_SI2_SI3_SI4_Chr4.txt &
	wait
	## SS
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o SS_Chr1-4_msmc_output ./multihetsep-Chr1/multihetsep_SS1_SS2_SS3_SS4_Chr1.txt ./multihetsep-Chr2/multihetsep_SS1_SS2_SS3_SS4_Chr2.txt ./multihetsep-Chr3/multihetsep_SS1_SS2_SS3_SS4_Chr3.txt ./multihetsep-Chr4/multihetsep_SS1_SS2_SS3_SS4_Chr4.txt

#### 5.2.7 Cross-coalesence for 16 haplotypes (2 populations á 4 individuals)
Need to perform 3 separate run + combineCrossCoal.py
1.    msmc2 -I 0,1,2,3,4,5,6,7 -o within1_msmc <input_chr1> <input_chr2> ...
2.    msmc2 -I 8,9,10,11,12,13,14,15 -o within2_msmc <input_chr1> <input_chr2> ...
3.    msmc2 -I 0-8,0-9,0-10,0-11,0-12,0-13,0-14,0-15, \
		1-8,1-9,1-10,1-11,1-12,1-13,1-14,1-15, \
		2-8,2-9,2-10,2-11,2-12,2-13,2-14,2-15, \
		3-8,3-9,3-10,3-11,3-12,3-13,3-14,3-15, \
		4-8,4-9,4-10,4-11,4-12,4-13,4-14,4-15, \
		5-8,5-9,5-10,5-11,5-12,5-13,5-14,5-15, \
		6-8,6-9,6-10,6-11,6-12,6-13,6-14,6-15, \
		7-8,7-9,7-10,7-11,7-12,7-13,7-14,7-15 -o across_msmc <input_chr1> <input_chr2> ...
4.    combineCrossCoal.py across_msmc.final.txt within1_msmc.final.txt within2_msmc.final.txt > combined12_msmc.final.txt

##### 1. Create multihetsep with 16 haplotypes

	#Cheops0
	module purge
	module load python/3.4.3

	cd /projects/ag-waldvogel/pophistory/CRIP/phasing/Chr1/phased

	# 4 inidividuals per population
	while read a b c d e f g h; do /home/lpettric/bin/msmc-tools//generate_multihetsep.py --mask=../$a"_Chr1_mask.bed.gz" \
                          --mask=../$b"_Chr1_mask.bed.gz" \
                          --mask=../$c"_Chr1_mask.bed.gz" \
                          --mask=../$d"_Chr1_mask.bed.gz" \
                          --mask=../$e"_Chr1_mask.bed.gz" \
                          --mask=../$f"_Chr1_mask.bed.gz" \
                          --mask=../$g"_Chr1_mask.bed.gz" \
                          --mask=../$h"_Chr1_mask.bed.gz" \
                          --mask=/projects/ag-waldvogel/pophistory/CRIP/masking/final-mask/mask_Chr1_145_50.bed.gz \
                          $a"_Chr1_phased_merged.vcf.gz" $b"_Chr1_phased_merged.vcf.gz" \
                          $c"_Chr1_phased_merged.vcf.gz" \
                          $d"_Chr1_phased_merged.vcf.gz" $e"_Chr1_phased_merged.vcf.gz" $f"_Chr1_phased_merged.vcf.gz" \
                          $g"_Chr1_phased_merged.vcf.gz" \
                          $h"_Chr1_phased_merged.vcf.gz"> /projects/ag-waldvogel/pophistory/CRIP/msmc2/multihetsep-Chr1/"multihetsep_"$a"-"$h"_joined_Chr1.txt"; done < /projects/ag-waldvogel/pophistory/CRIP/msmc2/list-populations-cc

For list-populations-cc see below

##### 2. Create within_msmc and across_msmc

list-populations-cc (original without numbers):
1.	MF1	MF2	MF3	MF4	MG2	MG3	MG4	MG5
2.	MF1	MF2	MF3	MF4	NMF1	NMF2	NMF3	NMF4
3.	MF1	MF2	MF3	MF4	SI1	SI2	SI3	SI4
4.	MF1	MF2	MF3	MF4	SS1	SS2	SS3	SS4
5.	MG2	MG3	MG4	MG5	NMF1	NMF2	NMF3	NMF4
6.	MG2	MG3	MG4	MG5	SI1	SI2	SI3	SI4
7.	MG2	MG3	MG4	MG5	SS1	SS2	SS3	SS4
8.	NMF1	NMF2	NMF3	NMF4	SI1	SI2	SI3	SI4
9.	NMF1	NMF2	NMF3	NMF4	SS1	SS2	SS3	SS4
10.	SI1	SI2	SI3	SI4	SS1	SS2	SS3	SS4

##### Example with MF1-MG5

	# Cheops1
	module purge

	cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/cross-coalescene

	# MF (within1) - MG (within2)
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0,1,2,3,4,5,6,7 -o ./run1/within1_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt &
	wait
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 8,9,10,11,12,13,14,15 -o ./run1/within2_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt &
	wait
	/home/lpettric/bin/msmc2/build/release/msmc2 -p 1*3+1*2+22*1+1*2+1*3 -t 5 -I 0-8,0-9,0-10,0-11,0-12,0-13,0-14,0-15,1-8,1-9,1-10,1-11,1-12,1-13,1-14,1-15,2-8,2-9,2-10,2-11,2-12,2-13,2-14,2-15,3-8,3-9,3-10,3-11,3-12,3-13,3-14,3-15,4-8,4-9,4-10,4-11,4-12,4-13,4-14,4-15,5-8,5-9,5-10,5-11,5-12,5-13,5-14,5-15,6-8,6-9,6-10,6-11,6-12,6-13,6-14,6-15,7-8,7-9,7-10,7-11,7-12,7-13,7-14,7-15 -o ./run1/across_msmc2_MF1-MG5 ../multihetsep-Chr1/multihetsep_MF1-MG5_joined_Chr1.txt ../multihetsep-Chr2/multihetsep_MF1-MG5_joined_Chr2.txt ../multihetsep-Chr3/multihetsep_MF1-MG5_joined_Chr3.txt ../multihetsep-Chr4/multihetsep_MF1-MG5_joined_Chr4.txt


##### 3. Create combined cross-coalescence

	#Cheops0
	module purge
	module load python/3.4.3 

	cd /projects/ag-waldvogel/pophistory/CRIP/msmc2/cross-coalescene/run1

	# MF-MG
	/home/lpettric/bin/msmc-tools/combineCrossCoal.py across_msmc2_MF1-MG5.final.txt within1_msmc2_MF1-MG5.final.txt within2_msmc2_MF1-MG5.final.txt > combinedMF1-MG5_msmc2.final.txt


##### => All results of combined and cross-coalescence in Excel, mean values per population calculated and plotted in R


#### 5.2.8 Account for uncertainities in coalescence





## 5.3 eSMC 

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

### 5.3.1 Create input file
eSMC2 analysis of CRIP
Create multihetsep per individual and chromosome

	#Cheops0
	module purge
	module load python/3.4.3

	cd /projects/ag-waldvogel/pophistory/CRIP/phasing/


	while read a b; do /home/lpettric/bin/msmc-tools/generate_multihetsep.py --mask=$b/$a"_"$b"_mask.bed.gz" \
                          --mask=/projects/ag-waldvogel/pophistory/CRIP/masking/final-mask/"mask_"$b"_145_50.bed.gz" \
                          $b/phased/$a"_"$b"_phased_merged.vcf.gz" > /projects/ag-waldvogel/pophistory/CRIP/esmc/multihetsep-files/"multihetsep_"$a"_"$b".txt"; done < /projects/ag-waldvogel/pophistory/CRIP/esmc/list-crip-chr
                          
                          

# ***Panagrolaimus kolymaensis***
# 1. First look at PKOL Genome

# 2. Pre-processing
