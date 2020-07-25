library(edgeR)
library(data.table)  
library(dplyr)

MasterFile <- read.csv("C:\\Users\\libin\\UCSF\\MECP2\\RNASeq\\TMM\\for_edgeR.csv")
sample_ID <- as.character(MasterFile$ID)
dir <- ("C:\\Users\\libin\\UCSF\\MECP2\\RNASeq\\TMM")
files <- file.path(dir, paste0(MasterFile$ID, "_rsem", ".genes.reformed.results"))

### Read in effective gene length of each libs
### 'length' is this transcript's sequence length (poly(A) tail is not counted). 
### 'effective_length' counts only the positions that can generate a valid fragment.
### If no poly(A) tail is added, 'effective_length' is equal to transcript length - mean fragment length + 1. 
### If one transcript's effective length is less than 1, this transcript's both effective length and abundance estimates are set to 0.
gene_length_all = lapply(files, read.delim2, stringsAsFactors=FALSE)
gene_length_all <- merge(gene_length_all[[1]], gene_length_all[[2]], by="gene_id")
row.names(gene_length_all) <- gene_length_all$gene_id
gene_length_all <- gene_length_all[c("effective_length.x", "effective_length.y")]
colnames(gene_length_all) <- sample_ID
scale_kb <- function(x) (x /1000)
# gene_length_all <- as.numeric(as.character(gene_length_all))
gene_length_all <- dplyr::mutate_all(gene_length_all, as.numeric)
gene_length_all_kb <- dplyr::mutate_all(gene_length_all, scale_kb)


DG_raw <- readDGE(files, columns=c(1,2), labels=sample_ID)
DG <- DGEList(counts = DG_raw[["counts"]])
# TMM normalization
DG_norm <- calcNormFactors(DG)

# rpkm_DG <- rpkm(DG_norm, gene.length = DG$genes$genes)

###
lib.size <- DG_norm$samples$lib.size
lib.size <- lib.size*DG_norm$samples$norm.factors
cpm_DG <- cpm.default(DG_norm,lib.size=lib.size,log=FALSE, prior.count=2)

RPKM_DG <- cpm_DG/gene_length_all_kb
row.names(RPKM_DG) <- rownames(cpm_DG)
write.table(RPKM_DG, file = "C:\\Users\\libin\\UCSF\\MECP2\\RNASeq\\TMM\\TMM_RPKM_DG", sep = "\t", quote = FALSE)
