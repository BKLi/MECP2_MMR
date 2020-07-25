#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=100G
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

export PATH=/wynton/scratch/bingkun/STAR:$PATH

ramGB=4

mkdir ${input}
cd "$input"
ln -s ../${input}_R1.trim_trimmed.fq.gz .
ln -s ../${input}_R2.trim_trimmed.fq.gz .
echo "$input"

STAR --genomeDir /wynton/scratch/bingkun/index_hg38_STAR \
    --readFilesIn "$input"_R1.trim_trimmed.fq.gz \
    --readFilesCommand zcat \
    --runThreadN 2 \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile COfile.txt \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate  \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM '${ramGB}'000000000
