#!/bin/bash
#$ -l h_rt=200:0:0
#$ -l mem_free=64G
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

export PATH=/shen/shenlabstore3/bingkun/miniconda3/bin:$PATH

juicertools=/wynton/home/shen/bingkun/Juicertools/juicer_tools_1.14.08.jar
input_hic=${input}.bsorted.pairs.hic
#CHR=19
res=5000
#output_matrix=${input}_${CHR}_${res}.matrix
maxmem=64G

for CHR in `seq 1 22` 
do 
    echo ${input}_"$CHR"_${res}.matrix
    java -Xmx$maxmem -Xms$maxmem -jar $juicertools dump observed NONE $input_hic $CHR $CHR BP $res ${input}_"$CHR"_${res}.matrix
done
