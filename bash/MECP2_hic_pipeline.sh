#!/bin/bash
#$ -l h_rt=60:0:0
#$ -l mem_free=120G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-3

i=$(($SGE_TASK_ID - 1))

#input_list=("DS001" "IJ18" \
#            "IJ14" "IJ39" "IJ255" \
#            "IJ12" "IJ13" "IJ243" \
#           "IJ19" "IJ20" \
#            "IJ239" "IJ142" "IJ251" \
#            "IJ240" "IJ144" "IJ253" \
#            "IJ173" "IJ143" "IJ250" \
#            "IJ174" "IJ145" "IJ252" "IJ254")

input_list=("IJ18" "IJ13" "IJ173")

input=${input_list[i]}

export PATH=/shen/shenlabstore3/bingkun/bingkun_home/bwa:$PATH
export PATH=/shen/shenlabstore3/bingkun/hg38_bowtie_bwa/bwa-index:$PATH

nThreads=2
index=/shen/shenlabstore3/bingkun/hg38_bowtie_bwa/bwa-index/hg38.fa
fastq1=${input}_R1.trimmed.fastq.gz
fastq2=${input}_R2.trimmed.fastq.gz
bwa mem -t $nThreads -SP5M $index $fastq1 $fastq2 | samtools view -Shb - > ${input}.bam
#!/bin/bash
#$ -l h_rt=60:0:0
#$ -l mem_free=60G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae
#$ -t 1-1

i=$(($SGE_TASK_ID - 1))

#input_list=("DS001" "IJ18" \
#           "IJ14" "IJ39" "IJ255" \
#          "IJ12" "IJ13" "IJ243" \
#            "IJ19" "IJ20" \
#           "IJ239" "IJ142" "IJ251" \
#            "IJ240" "IJ144" "IJ253" \
#            "IJ143" "IJ250" \
#            "IJ174" "IJ145" "IJ252" "IJ254")

input_list=("IJ173")

input=${input_list[i]}

mkdir ${input}
cd ${input}
ln -s ../${input}.bam .

samtools sort -o ${input}.sorted.bam -T ${input}.tmp.bam ${input}.bam
samtools index ${input}.sorted.bam
#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=60G
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
            "IJ12" "IJ13" "IJ243" \
            "IJ19" "IJ20" \
            "IJ239" "IJ251" \
            "IJ144" "IJ253" \
            "IJ143" "IJ250" "IJ173" \
            "IJ174" "IJ145" "IJ252" "IJ254")

input=${input_list[i]}

source ~/.bashrc
conda activate 4DN

# N_THREADS=8
PAIRSAM=${input}.sam.pairs.gz
MARKED_PAIRSAM=${input}.marked.sam.pairs.gz
LOSSLESS_BAM=${input}.lossless.bam
UNMAPPED_SAMPAIRS=${input}.unmapped.sam.pairs.gz
TEMPFILE=${input}.temp.gz
TEMPFILE1=${input}.temp1.gz
DEDUP_PAIRS=${input}.dedup.pairs.gz
CHR_SIZES=/shen/shenlabstore3/bingkun/hg38_for_hicpro_hicup/hg38.chrom.sizes


echo "Parsing & Sorting ..."
echo input:${input}.sorted.bam
samtools view -h ${input}.sorted.bam | pairtools parse -c ${CHR_SIZES} --add-columns mapq | pairtools sort --compress-program gzip --output ${PAIRSAM}
echo output:${PAIRSAM}
echo "=== parsing & sorting finished! ==="


echo "Marking and removing PCR duplicates..."
echo input:${PAIRSAM}
pairtools dedup --mark-dups --output-dups - --output-unmapped - --output ${MARKED_PAIRSAM} ${PAIRSAM}
pairix ${MARKED_PAIRSAM} # sanity check
echo output:${MARKED_PAIRSAM}
echo "=== Deduplication finished! ==="


echo "Generating lossless bam..."
echo input:${MARKED_PAIRSAM}
pairtools split --output-sam ${LOSSLESS_BAM} ${MARKED_PAIRSAM}
echo output:${LOSSLESS_BAM}
echo "=== Lossless bam generated! ==="


echo "Filtering ... Selecting UU, UR, RU reads"
echo input:$MARKED_PAIRSAM
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
--output-rest ${UNMAPPED_SAMPAIRS} \
--output ${TEMPFILE} \
$MARKED_PAIRSAM
echo output:${UNMAPPED_SAMPAIRS} ${TEMPFILE}
echo "=== Filtering done! ==="


echo "Spliting filtered file..."
echo input:${TEMPFILE}
pairtools split --output-pairs ${TEMPFILE1} ${TEMPFILE}
echo output:${TEMPFILE1}
echo "=== Splitting done! ==="


echo "Generating final pairs file..."
echo input:${TEMPFILE1}
pairtools select 'True' --chrom-subset ${CHR_SIZES} -o ${DEDUP_PAIRS} ${TEMPFILE1}
pairix ${DEDUP_PAIRS}
echo output:${DEDUP_PAIRS}
echo "=== Final pairs file generated! ==="


echo "removing temp files..."
rm ${TEMPFILE1}
rm ${TEMPFILE}
rm ${MARKED_PAIRSAM}
rm ${PAIRSAM}

echo "=== ALL DONE! ==="
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
input_pairs=${input}.dedup.pairs.gz
output_hic=${input}.dedup.pairs.hic
chromsizefile=/shen/shenlabstore3/bingkun/hg38_for_hicpro_hicup/hg38.chrom.sizes
mapqfilter=30
maxmem=64G

java -Xmx$maxmem -Xms$maxmem -jar $juicertools pre -n $input_pairs $output_hic $chromsizefile 
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
input_hic=${input}.dedup.ff.pairs.hic
CHR=19
res=5000
output_matrix=${input}_${CHR}_${res}.matrix
maxmem=64G

java -Xmx$maxmem -Xms$maxmem -jar $juicertools dump observed NONE $input_hic $CHR $CHR BP $res $output_matrix
