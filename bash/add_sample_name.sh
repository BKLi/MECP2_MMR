#!/bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -t 1-35

i=$(($SGE_TASK_ID - 1))

export PATH=/shen/shenlabstore3/bingkun/bingkun_home/FastQC:$PATH

input_list=("LM036" "LM037" \
            "LM044" "LM045" \
            "LM065" "LM066" "LM067" "LM068" "LM069" \
            "LM070" "LM071" "LM072" "LM073" "LM074" "LM075" "LM076" \
            "LM118" "LM119" \
            "LM120" "LM121" "LM122" \
            "LM124" "LM125" "LM126" "LM134" \
            "LMA017" "LMA035" \
            "LMA045" "LMA046" "LMA047" "LMA048" "LMA049" \
            "LMA050" "LMA051" "LMA052")

input=${input_list[i]}

cd ${input}
sample_name=$(echo "${PWD##*/}")
echo $sample_name

file=Aligned.toTranscriptome.out.bam

mv $file "$sample_name"_$file
