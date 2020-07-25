#!/bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-17

i=$(($SGE_TASK_ID - 1))

export PATH=/shen/shenlabstore3/bingkun/TrimGalore/TrimGalore-0.6.0:$PATH

input_list=('LM065' 'LM069' 'LM066' 'LM071' \
            'LM067' 'LM045' 'LM068' 'LMA017' \
            'LM075' 'LM036' 'LM073' 'LM044' \
            'LM037' 'LM074' 'LM070' 'LM076' 'LM072')

input=${input_list[i]}
trim_galore -q 20 --length 20 --stringency 3 --trim-n ${input}_R1.trim.fastq.gz

