#!/bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-18

i=$(($SGE_TASK_ID - 1))

export PATH=/shen/shenlabstore3/bingkun/TrimGalore/TrimGalore-0.6.0:$PATH

input_list=('LM118' 'LM119' 'LM120' \
            'LM121' 'LM122' 'LM124' \
            'LM125' 'LM126' 'LM134' \
            'LMA035' 'LMA045' 'LMA046' \
            'LMA047' 'LMA048' 'LMA049' \
            'LMA050' 'LMA051' 'LMA052')

input=${input_list[i]}
trim_galore --paired -q 20 --length 20 --stringency 3 --trim-n ${input}_R1.trim.fastq.gz ${input}_R2.trim.fastq.gz

