#!/bin/bash
#$ -l h_rt=48:0:0
#$ -l mem_free=20G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
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

export PATH=/shen/shenlabstore3/bingkun/bingkun_home/bwa:$PATH
export PATH=/shen/shenlabstore3/bingkun/miniconda3/envs/hicup/bin:$PATH
export PATH=/shen/shenlabstore3/bingkun/miniconda3/bin:$PATH

# samtools sort -T ${input}_tmp -o ${input}.sorted.bam ${input}_Aligned.toTranscriptome.out.bam
# samtools index ${input}.sorted.bam
bamCoverage --bam ${input}.sorted.bam -o ${input}.bw
