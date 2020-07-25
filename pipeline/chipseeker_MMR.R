library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(UpSetR)
library(ggplot2)
library(grid)

hg38_txdb <- makeTxDbFromGFF("C:\\Users\\libin\\UCSF\\hg38_general\\gencode.v32.annotation.gtf", format = c("gtf"))
MMR_df <- read.table("C:\\Users\\libin\\UCSF\\MMR\\mmr_elements.txt", sep = "\t", header = TRUE)
MMR_df$ID <- paste("element",seq.int(nrow(MMR_df)), sep="_")
write.table(MMR_df, file="C:\\Users\\libin\\UCSF\\MMR\\mmr_elements_withID", sep="\t", quote = FALSE, row.names = FALSE)
MMR_grange <- makeGRangesFromDataFrame(MMR_df, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = "chr", start.field = "Start", end.field = "End" )

MMR_ann_EXON <- annotatePeak(MMR_grange, tssRegion = c(-500,500), TxDb = hg38_txdb, 
                             genomicAnnotationPriority = c("Exon", "Promoter", "Intron", "Intergenic", "5UTR", "3UTR",  "Downstream"))
MMR_ann_df_EXON <- as.data.frame(MMR_ann_EXON)
write.table(MMR_ann_df_EXON, file=sprintf("C:\\Users\\libin\\UCSF\\MMR\\%s.MMR.chipseek.output.bed", "exon"), quote = FALSE,sep = "\t", row.names = FALSE)
###
MMR_ann_INTRON <- annotatePeak(MMR_grange, tssRegion = c(-500,500), TxDb = hg38_txdb, 
                             genomicAnnotationPriority = c("Intron", "Exon", "Promoter", "5UTR", "3UTR", "Intergenic", "Downstream"))
MMR_ann_df_INTRON <- as.data.frame(MMR_ann_INTRON)
write.table(MMR_ann_df_INTRON, file=sprintf("C:\\Users\\libin\\UCSF\\MMR\\%s.MMR.chipseek.output.bed", "intron"), quote = FALSE,sep = "\t", row.names = FALSE)
###
MMR_ann_promoter <- annotatePeak(MMR_grange, tssRegion = c(-500,500), TxDb = hg38_txdb, 
                               genomicAnnotationPriority = c("Promoter", "Intron", "Exon", "5UTR", "3UTR", "Intergenic", "Downstream"))
MMR_ann_df_promoter <- as.data.frame(MMR_ann_promoter)
write.table(MMR_ann_df_promoter, file=sprintf("C:\\Users\\libin\\UCSF\\MMR\\%s.MMR.chipseek.output.bed", "promoter"), quote = FALSE,sep = "\t", row.names = FALSE)
###
MMR_ann_intergenic <- annotatePeak(MMR_grange, tssRegion = c(-500,500), TxDb = hg38_txdb, 
                                 genomicAnnotationPriority = c("Intergenic", "Promoter", "Intron", "Exon", "5UTR", "3UTR", "Downstream"))
MMR_ann_df_intergenic <- as.data.frame(MMR_ann_intergenic)
write.table(MMR_ann_df_intergenic, file=sprintf("C:\\Users\\libin\\UCSF\\MMR\\%s.MMR.chipseek.output.bed", "intergenic"), quote = FALSE,sep = "\t", row.names = FALSE)
