#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=60G
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
            "IJ143" "IJ250" "IJ173" \
            "IJ174" "IJ145" "IJ252" "IJ254")

# input_list=("IJ18" "IJ13" "IJ173")

input=${input_list[i]}

source ~/.bashrc
conda activate 4DN

# N_THREADS=8
PAIRSAM=${input}.sam.pairs.gz
MARKED_PAIRSAM=${input}.marked.sam.pairs.gz
LOSSLESS_BAM=${input}.lossless.bam
UNMAPPED_SAMPAIRS=${input}.unmapped.sam.pairs.gz
TEMPFILE=temp.gz
TEMPFILE1=temp1.gz
DEDUP_PAIRS=${input}.dedup.pairs.gz
CHR_SIZES=/shen/shenlabstore3/bingkun/hg38_for_hicpro_hicup/hg38.chrom.sizes

# PARSE & SORT
samtools view -h ${input}.sorted.bam | pairtools parse -c ${CHR_SIZES} --add-columns mapq | pairtools sort --compress-program gzip --output ${PAIRSAM}

# REMOVE PCR DUPLICATES
pairtools dedup --mark-dups --output-dups - --output-unmapped - --output ${MARKED_PAIRSAM} ${PAIRSAM}
pairix ${MARKED_PAIRSAM} # sanity check

# generate lossless bam
pairtools split --output-sam ${LOSSLESS_BAM} ${MARKED_PAIRSAM}

# FILTERING
# Select UU, UR, RU reads
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
--output-rest ${UNMAPPED_SAMPAIRS} \
--output ${TEMPFILE} \
$MARKED_PAIRSAM

pairtools split --output-pairs ${TEMPFILE1} ${TEMPFILE}

pairtools select 'True' --chrom-subset ${CHR_SIZES} -o ${DEDUP_PAIRS} ${TEMPFILE1}

pairix ${DEDUP_PAIRS}
