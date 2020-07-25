#!/bin/bash
#$ -l h_rt=60:0:0
#$ -l mem_free=120G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-3

i=$(($SGE_TASK_ID - 1))

#input_list=("DS001" "IJ18" \
#            "IJ14" "IJ39" "IJ255" \
#            "IJ12" "IJ13" "IJ243" \
#           "IJ19" "IJ20" \
#            "IJ239" "IJ142" "IJ251" \
#            "IJ240" "IJ144" "IJ253" \
#            "IJ173" "IJ143" "IJ250" \
#            "IJ174" "IJ145" "IJ252" "IJ254")

input_list=("IJ18" "IJ13" "IJ173")

input=${input_list[i]}

export PATH=/shen/shenlabstore3/bingkun/bingkun_home/bwa:$PATH
export PATH=/shen/shenlabstore3/bingkun/hg38_bowtie_bwa/bwa-index:$PATH

nThreads=2
index=/shen/shenlabstore3/bingkun/hg38_bowtie_bwa/bwa-index/hg38.fa
fastq1=${input}_R1.trimmed.fastq.gz
fastq2=${input}_R2.trimmed.fastq.gz
bwa mem -t $nThreads -SP5M $index $fastq1 $fastq2 | samtools view -Shb - > ${input}.bam
