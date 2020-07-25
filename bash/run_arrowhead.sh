#!/bin/bash
#$ -l h_rt=200:0:0
#$ -l mem_free=64G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-21

i=$(($SGE_TASK_ID - 1))

input_list=("DS001" "IJ18" \
            "IJ14" "IJ39" "IJ255" \
            "IJ12" "IJ243" \
            "IJ19" "IJ20" \
            "IJ239" "IJ142" "IJ251" \
            "IJ240" "IJ144" "IJ253" \
            "IJ143" "IJ250"  \
            "IJ174" "IJ145" "IJ252" "IJ254")
# input_list=("IJ13" "IJ173")
input=${input_list[i]}

export PATH=/shen/shenlabstore3/bingkun/miniconda3/bin:$PATH

juicertools=/wynton/home/shen/bingkun/Juicertools/juicer_tools_1.14.08.jar
input_hic=${input}.bsorted.pairs.hic
res=50000
output_domain=${input}_${res}.2D
maxmem=64G
thr=8

java -Xmx$maxmem -Xms$maxmem -jar $juicertools arrowhead -r $res --threads $thr -k KR $input_hic $output_domain
