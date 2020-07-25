#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=80G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m a
#$ -t 1-35

i=$(($SGE_TASK_ID - 1))

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

cd "$input"
echo "$input"

export PATH=/shen/shenlabstore3/bingkun/RSEM/RSEM-build/bin:$PATH

# use toTranscriptome for RSEM quantification
rsem-sam-validator Aligned.toTranscriptome.out.bam

