setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(DESeq2)
library(pheatmap)

# Obtain dds from DESeq2.R
# Gene orthologs from TAIR
# GO annotation from entap
GO_anno <- read.csv("../results/annotations_GO.tsv", sep = "\t", header = T)
GO_anno$Query.Sequence <- gsub("ID=","",GO_anno$Query.Sequence)
GO_anno <- GO_anno[order(GO_anno$Query.Sequence),]
# get the transformed gene matrix
vsd <- vst(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("Condition","Stage")])
names(df) <- c("Stage","Condition")
vsd_ordered <- vsd[,c(paste0("S",1:6), "S7","S12", paste0("S",8:11),paste0("S",13:18))]

# photosynthesis genes
photosyn_blast <- read.table("../results/Vd_photosyn.tsv", sep = "\t", header = F)
length(unique(photosyn_blast$V1)) # count how many genes - 1,218 genes in total
photosyn_blast$V1 <- gsub("ID=","",photosyn_blast$V1) # remove "ID="
KEGG_photo <- GO_anno[grep("196|195",GO_anno$EggNOG.KEGG.Terms),]
KEGG_photo <- KEGG_photo[-which(is.na(KEGG_photo$EggNOG.Description)),]
# filter blast results
photosyn_blast_filter <- photosyn_blast[which(photosyn_blast$V3 > 30),]
# extract DEGs
GvR_sig <- read.csv("../results/DE_GvsR.csv", header = T)
select <- rownames(vsd) %in% intersect(KEGG_photo$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))
df$Condition <- factor(df$Condition, levels = c("Green","Early_pink","Late_pink","Ripening"))
ann_colors = list(
  Stage = c(Green="#009E73", Early_pink="#E69F00",Late_pink="#CC79A7",Ripening="#0072B2")
)
data <- assay(vsd_ordered)[select,]
rownames(data) <- KEGG_photo$EggNOG.Description[KEGG_photo$Query.Sequence %in% intersect(KEGG_photo$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))]
pheatmap(data, filename = "../img/Photosynthesis_heatmap.pdf", width = 10, height = 4.5,
         cluster_rows=T, show_rownames=T,show_colnames = F,
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])

# flavonoid genes
flavo_blast <- read.table("../results/Vd_flavo.tsv", sep = "\t", header = F)
flavo_blast$V1 <- gsub("ID=","",flavo_blast$V1) # remove "ID="
flavo_blast_filter <- flavo_blast[which(flavo_blast$V3 > 60),]
KEGG_antho <- GO_anno[grep("942|941",GO_anno$EggNOG.KEGG.Terms),]
#KEGG_antho <- KEGG_antho[-which(is.na(KEGG_antho$EggNOG.Description)),]
# extract DEGs
GvEP_sig <- read.csv("../results/DE_GvsEP.csv", header = T)
GvLP_sig <- read.csv("../results/DE_GvsLP.csv", header = T)
select <- rownames(vsd) %in% intersect(KEGG_antho$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))
data <- assay(vsd_ordered)[select,]
rownames(data) <- KEGG_antho$EggNOG.Description[KEGG_antho$Query.Sequence %in% intersect(KEGG_antho$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))]
rownames(data)[grep("primary product",rownames(data))] <- "chalcone synthase"
rownames(data)[19] <- "Methylates caffeoyl CoA to feruloyl CoA"
pheatmap(data, filename = "../img/Flavo-antho_heatmap_scaled.pdf", width = 10, height = 5,
         cluster_rows=T, show_rownames=T,show_colnames = F,scale="row",
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])

# ROS genes
GO_ros <- GO_anno[grep("GO:0006979",GO_anno$EggNOG.GO.Biological),]
GO_ros <- GO_ros[-which(is.na(GO_ros$EggNOG.Description)),]
# extract DEGs
LPvR_sig <- read.csv("../results/DE_LvsR.csv", header = T)
EPvR_sig <- read.csv("../results/DE_EPvsR.csv", header = T)
select <- rownames(vsd) %in% intersect(GO_ros$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))
data <- assay(vsd_ordered)[select,]
rownames(data) <- GO_ros$EggNOG.Description[GO_ros$Query.Sequence %in% intersect(GO_ros$Query.Sequence, c(GvR_sig$X,GvLP_sig$X,LPvR_sig$X,GvEP_sig$X))]
pheatmap(data, filename = "../img/ROS_heatmap.pdf", width = 5, height = 5,
         cluster_rows=T, show_rownames=F,show_colnames = F,
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])
