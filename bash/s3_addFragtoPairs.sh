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

input=${input_list[i]}

source ~/.bashrc
conda activate 4DN

input_pair=${input}.dedup.pairs.gz 
output_pair=${input}.dedup.ff.pairs
restriction_file=/wynton/scratch/bingkun/MECP2/HiC/hg38_DpnII.txt

which fragment_4dnpairs.pl

zcat $input_pair | fragment_4dnpairs.pl -a - $output_pair $restriction_file
bgzip -f $output_pair
pairix -f "$output_pair".gz
