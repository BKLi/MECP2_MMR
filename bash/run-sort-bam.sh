#!/bin/bash
#$ -l h_rt=60:0:0
#$ -l mem_free=60G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-22

i=$(($SGE_TASK_ID - 1))

input_list=("DS001" "IJ18" \
            "IJ14" "IJ39" "IJ255" \
            "IJ12" "IJ13" "IJ243" \
            "IJ19" "IJ20" \
            "IJ239" "IJ142" "IJ251" \
            "IJ240" "IJ144" "IJ253" \
            "IJ143" "IJ250" \
            "IJ174" "IJ145" "IJ252" "IJ254")

# input_list=("IJ18" "IJ13" "IJ173")

input=${input_list[i]}

mkdir ${input}
cd ${input}
ln -s ../${input}.bam .

samtools sort -o ${input}.sorted.bam -T ${input}.tmp.bam ${input}.bam
samtools index ${input}.sorted.bam
