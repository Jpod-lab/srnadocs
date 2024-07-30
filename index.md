## Samples
initial fastq files smallRNA


10_12d_4hrA.fastq.gz 17_12d_24hrA.fastq.gz 36_4d_4hrA.fastq.gz
5_12d_t0.fastq.gz 73_20d_2hrR.fastq.gz 101_D2_t0.fastq.gz
18_12d_24hrA.fastq.gz 37_4d_4hrA.fastq.gz 52_20d_t0.fastq.gz
74_20d_2hrR.fastq.gz 102_D2_t0.fastq.gz 19_12d_2hrR.fastq.gz
38_4d_4hrA.fastq.gz 53_20d_t0.fastq.gz 76_20d_24hrR.fastq.gz
104_D2_4hrA.fastq.gz 21_12d_2hrR.fastq.gz 39_4d_24hrA.fastq.gz
55_20d_t0.fastq.gz 77_20d_24hrR.fastq.gz 108_D2_24hrA.fastq.gz
2_12d_t0.fastq.gz 40_4d_24hrA.fastq.gz 57_20d_2hrA.fastq.gz
79_20d_24hrR.fastq.gz 110_D2_2hrR.fastq.gz 22_12d_2hrR.fastq.gz
41_4d_24hrA.fastq.gz 59_20d_2hrA.fastq.gz 80_20d_24hrR.fastq.gz
11_12d_4hrA.fastq.gz 23_12d_2hrR.fastq.gz 42_4d_24hrA.fastq.gz
60_20d_2hrA.fastq.gz 81_D2_t0.fastq.gz 111_D2_2hrR.fastq.gz
25_12d_24hrR.fastq.gz 43_4d_2hrR.fastq.gz 61_20d_2hrA.fastq.gz
82_D2_4hrA.fastq.gz 112_D2_24hrR.fastq.gz 28_12d_24hrR.fastq.gz
44_4d_2hrR.fastq.gz 6_12d_t0.fastq.gz 86_D2_t0.fastq.gz
1_12d_t0.fastq.gz 29_12d_24hrR.fastq.gz 45_4d_2hrR.fastq.gz
63_20d_6hrA.fastq.gz 89_D2_4hrA.fastq.gz 114_D2_24hrA.fastq.gz
30_12d_24hrR.fastq.gz 46_4d_2hrR.fastq.gz 64_20d_6hrA.fastq.gz
91_D2_4hrA.fastq.gz 118_D2_24hrR.fastq.gz 31_4d_t0.fastq.gz
47_4d_24hrR.fastq.gz 65_20d_6hrA.fastq.gz 92_D2_24hrA.fastq.gz
119_D2_24hrR.fastq.gz 32_4d_t0.fastq.gz 48_4d_24hrR.fastq.gz
67_20d_6hrA.fastq.gz 94_D2_24hrA.fastq.gz 12_12d_4hrA.fastq.gz
33_4d_t0.fastq.gz 49_4d_24hrR.fastq.gz 70_20d_2hrR.fastq.gz
96_D2_2hrR.fastq.gz 13_12d_24hrA.fastq.gz 34_4d_t0.fastq.gz
50_4d_24hrR.fastq.gz 7_12d_4hrA.fastq.gz 97_D2_2hrR.fastq.gz
16_12d_24hrA.fastq.gz 35_4d_4hrA.fastq.gz 51_20d_t0.fastq.gz
72_20d_2hrR.fastq.gz 98_D2_24hrR.fastq.gz

## Other input files
1. File for adaptor contamination (.fa/.txt)
2. Reference genome with mito (GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.fna)

## Install and load software

module load R/4.4.0/gcc-12.1.0 
module load gcc-11.1.0 
conda create -n fastqc fastqc 
conda activate fastqc 
conda install -n multiqc

## Terminal multiplexer

tmux new -s fastqc

## FASTQC 

fastqcr.R

if (!require(\"fastqcr\")) { install.packages(\"fastqcr\",
repos=\"http://cran.us.r-project.org\") library(fastqcr) } else
{library(fastqcr)} fastqc(fq.dir =
\"/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/copied_old_infiles/init_fastq_files\",
qc.dir=
\"/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/fastqc/init_fastqc\",
threads=20)

Rscript fastqcr.R

## Trimming 

module load Biosciences/Trimmomatic/0.39 Old script -threads 6 -phred33
-trimlog
/disk/scratch/Podrab/Claire/TrimRes_alim_anox_092716/10_12d_4hrA_trim.txt
/disk/scratch/Podrab/Claire/alim_anox_concat_092716/10_12d_4hrA.fastq.gz
/disk/scratch/Podrab/Claire/TrimRes_alim_anox_092716/10_12d_4hrA_trimm.fastq
ILLUMINACLIP:/vol/share/podrabsj_lab/Amielynn/adapter_contamination_sequences_AR3.txt:2:30:5:1:true
SLIDINGWINDOW:5:15 LEADING:20 TRAILING:20 MINLEN:15

FILES=\"/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/copied_old_infiles/init_fastq_files/\*.fastq.gz\"
for f in \$FILES do b=\`basename \$f\` echo Running trimming for the
file \$b c=\${b::-9} o=\"\$c.trim.fastq\" c=\"\$c.log.txt\" java -jar
/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/trim/Trimmomatic-0.39/trimmomatic-0.39.jar
SE -threads 20 -phred33 -trimlog \$c \$f \$o
/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/copied_old_infiles/adapter_contamination_sequences_AR3.txt:2:30:5:1:true
SLIDINGWINDOW:5:15 LEADING:20 TRAILING:20 MINLEN:15 done #done

tmux new -s fastqc tmux a -t fastqc

## FASTQC 

tmux new -s fastqc tmux a -t fastqc

fastqcr.R

if (!require(\"fastqcr\")) { install.packages(\"fastqcr\",
repos=\"http://cran.us.r-project.org\") library(fastqcr) } else
{library(fastqcr)} fastqc(fq.dir =
\"/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/trim\",
qc.dir=\"/disk/bioscratch/Podrab_lab/gazal/sRNA_gazal/fastqc/trim_fastqc\",
threads=20)

## Contributors
Gazal Kalyan \
Amie Romney \
Claire Riggs \
Jason Podrabsky 
