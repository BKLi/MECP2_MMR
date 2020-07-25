#!/bin/bash
#$ -l h_rt=100:0:0
#$ -l mem_free=100G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -M libingkun1997@gmail.com
#$ -m ae

unset PYTHONPATH
export PYTHONUSERSITE=1
export TMPDIR=/wynton/scratch/bingkun/tmp

export PATH=/wynton/scratch/bingkun/.bds:$PATH
export PATH=/wynton/scratch/bingkun/hg38_atac:$PATH
export PATH=/shen/shenlabstore3/bingkun/atac_dnase_pipelines:$PATH
export PATH=/shen/shenlabstore3/bingkun/miniconda3/envs/bds_atac/extra/phantompeakqualtools:$PATH 

export PATH=/shen/shenlabstore3/bingkun/miniconda3/envs/bds_atac/bin:$PATH
export PATH=/shen/shenlabstore3/bingkun/miniconda3/envs/bds_atac_py3/bin:$PATH

/wynton/scratch/bingkun/.bds/bds /shen/shenlabstore3/bingkun/atac_dnase_pipelines/atac.bds \
-type atac-seq -species hg38 -auto_detect_adapter -se -sge \
-nth 6 -mapq_thresh 30 -memory 40G \
-mem_bwt2 60G -mem_trim 30G -mem_dedup 60G -mem_shuf 20G -mem_macs2 60G -mem_ataqc 30G \
-enable_idr -ENCODE3 \
-fastq1 LM062_R1.trim.fastq.gz \
-fastq2 LM085_R1.trim.fastq.gz \
-fastq3 LM064_R1.trim.fastq.gz \
-out_dir atac_0_2E10 \
-chrsz /wynton/scratch/bingkun/hg38_atac/hg38/hg38.chrom.sizes \
-seq /wynton/scratch/bingkun/hg38_atac/hg38/seq \
-gensz hs \
-bwt2_idx /wynton/scratch/bingkun/hg38_atac/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
-ref_fa /wynton/scratch/bingkun/hg38_atac/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
-blacklist /wynton/scratch/bingkun/hg38_atac/hg38/hg38.blacklist.bed.gz \
-tss_enrich /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/hg38_gencode_tss_unique.bed.gz \
-dnase /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz \
-prom /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz \
-enh /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz \
-reg2map /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz \
-reg2map_bed /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz \
-roadmap_meta /wynton/scratch/bingkun/hg38_atac/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt





