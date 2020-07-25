library(DESeq2)
library(tximport)
library(pheatmap)
library(RColorBrewer)

MasterFile <- read.csv("C:\\Users\\libin\\UCSF\\MECP2\\for_deseq_0115.csv")
dir <- ("C:\\Users\\libin\\UCSF\\MECP2\\RNASeq")
files <- file.path(dir, paste0(MasterFile$ID, "_rsem", ".genes.results"))
names(files) <- MasterFile$ID
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length== 0] <- 1
ddsTxi <- DESeqDataSetFromTximport(txi.rsem, MasterFile, ~condition)
rsemDESeq <- DESeq(ddsTxi)
rsem_res <- results(rsemDESeq)

tsf.rsem <- vst(rsemDESeq)
# fpkm_rsem <- fpkm(rsemDESeq)
# rsem_res = results(rsemDESeq)
gene_df <- as.data.frame(colData(rsemDESeq)[,"condition"])
# rename column names
colnames(gene_df) <- "cell type"
#
# select_genes <- order(rowMeans(counts(rsemDESeq, normalized=TRUE)),decreasing = TRUE)[1:50]
# cortical_genes <- order((rowMeans(subset(counts(rsemDESeq, normalized=TRUE),select=c("LM007", "LM008"))) - rowMeans(subset(counts(rsemDESeq, normalized=TRUE),select=c("LM014", "LM015", "LMA024", "LMA025", "LMA028", "LMA033", "LMA029", "LMA032", "LMA063")))),decreasing = T)[1:20]
# hippo_sig_genes <- c("PROX1","TBR1",'MAP2',"RBFOX3", "TUBB3", "OCT4","PAX6","NEUROD1","DCX", "FOXG1")
# hippo_ensembl <- queryMany(hippo_sig_genes, scopes = "symbol", fields='ensembl.gene', "species"="human")
# hippo_ensembl_ID <- hippo_ensembl$ensembl.gene
# hippo_ID_pos <- match(hippo_ensembl_ID, rownames(tsf.rsem))
# hippo_df <- subset(assay(tsf.rsem), select=c("LMA029","LMA032","LMA063"))
# motor_df_geneName <- motor_df[motor_ID_pos,]
# rownames(cor_df_geneName) <- cortical_ensembl$query
# pheatmap(assay(tsf.rsem)[hippo_ID_pos,], cluster_cols = TRUE, show_rownames = TRUE, annotation_col = gene_df) 

# pheatmap(assay(tsf.rsem)[select_genes,], cluster_cols = TRUE, show_rownames = TRUE, annotation_col = gene_df) 

sampledDists <- dist(t(assay(tsf.rsem)), method = "euclidean")
sampleDistMatrix <- as.matrix(sampledDists)
rownames(sampleDistMatrix) <- MasterFile$ID
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampledDists, clustering_distance_cols=sampledDists, col=colors)
