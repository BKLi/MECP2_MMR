#!/bin/bash
#$ -l h_rt=10:0:0
#$ -l mem_free=40G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-23

i=$(($SGE_TASK_ID - 1))

input_list=("DS001" "IJ18" \
            "IJ14" "IJ39" "IJ255" \
            "IJ12" "IJ13" "IJ243" \
            "IJ19" "IJ20" \
            "IJ239" "IJ142" "IJ251" \
            "IJ240" "IJ144" "IJ253" \
            "IJ173" "IJ143" "IJ250" \
            "IJ174" "IJ145" "IJ252" "IJ254")

input=${input_list[i]}

export PATH=/wynton/scratch/bingkun/FastQC:$PATH

# fastqc ${input}_R1.trimmed.fastq.gz
fastqc ${input}_R2.trimmed.fastq.gz
