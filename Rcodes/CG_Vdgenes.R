setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

CG_Genes <- read.table("../results/CG_Vc2VdGenes.tsv", header = F, sep = "\t")
Vd_CG <- unique(CG_Genes$V2)
Vd_CG_annot <- data.frame(ID = Vd_CG, label=c("CYP791A-like",rep("CYP71E1-like",2),"UGT85B1-like",
                                              "dhurrinease-like",rep("hydroxynitrile lyase-like",2),
                                              "NIT4-like"))
Vd_CG_annot$ID <- gsub("ID=","",Vd_CG_annot$ID)

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
select <- rownames(vsd) %in% Vd_CG_annot$ID
data <- assay(vsd_ordered)[select,]
rownames(data) <- paste0(Vd_CG_annot[order(Vd_CG_annot$ID),]$label,":VaDar_",Vd_CG_annot[order(Vd_CG_annot$ID),]$ID)
pheatmap(data, filename = "../img/Vd_CG_heatmap_scaled.pdf", width = 10, height = 3,
         cluster_rows=T, show_rownames=T,show_colnames = F,scale = "row",
         cluster_cols=F, annotation_col=df[,-2, drop=F],border_color=F, annotation_colors = ann_colors[1])
