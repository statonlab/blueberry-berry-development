setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(DESeq2)
library(pheatmap)

# Vcorymbosum anthocyanin genes
antho_Genes <- read.csv("../dataInfo/Vc_anthocy_genes.csv", header = T)
RBH <- read.table("../results/reciprocal_best_hits.txt", header = F, sep = "\t")
antho_TF <- read.csv("../dataInfo/Vc_anthocy_TF.csv", header = T)

head(RBH)
head(antho_Genes)
antho_Genes$Gene.Name <- gsub("maker-|augustus_masked-|snap_masked-", "", antho_Genes$Gene.Name)
antho_TF$Genes <- gsub("maker-|augustus_masked-|snap_masked-", "", antho_TF$Genes)
RBH$V1 <- gsub("-mRNA-1-mRNA-1","",RBH$V1)

antho_RBH <- RBH[RBH$V1 %in% antho_Genes$Gene.Name,]
antho_RBH$V2 <- gsub("ID=","",antho_RBH$V2)
antho_RBH <- merge(antho_RBH, antho_Genes, by.x = "V1", by.y = "Gene.Name", all.x=T)
antho_TF_RBH <- RBH[RBH$V1 %in% antho_TF$Genes,]
antho_TF_RBH <- merge(antho_TF_RBH, antho_TF, by.x = "V1", by.y = "Genes", all.x=T)
antho_TF_RBH$V2 <- gsub("ID=","",antho_TF_RBH$V2)

# Obtain dds from DESeq2.R
# get the transformed gene matrix
ann_colors = list(
  Stage = c(Green="#009E73", Early_pink="#E69F00",Late_pink="#CC79A7",Ripening="#0072B2")
)
vsd <- vst(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("Condition","Stage")])
names(df) <- c("Stage","Condition")
vsd_ordered <- vsd[,c(paste0("S",1:6), "S7","S12", paste0("S",8:11),paste0("S",13:18))]
# plot anthocyanin biosynthesis genes
select <- rownames(vsd) %in% antho_RBH$V2
data <- assay(vsd_ordered)[select,]
rownames(data) <- antho_RBH[order(antho_RBH$V2),]$Homolog
pheatmap(data, filename = "../img/VcRbh_anthoGenes_heatmap.pdf", width = 10, height = 4,
         cluster_rows=T, show_rownames=T,show_colnames = F,
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])

# plot anthocyanin TF genes
select <- rownames(vsd) %in% antho_TF_RBH$V2
data <- assay(vsd_ordered)[select,]
rownames(data) <- antho_TF_RBH[order(antho_TF_RBH$V2),]$label
pheatmap(data, filename = "../img/VcRbh_anthoTF_heatmap.pdf", width = 10, height = 3,
         cluster_rows=T, show_rownames=T,show_colnames = F,
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])
