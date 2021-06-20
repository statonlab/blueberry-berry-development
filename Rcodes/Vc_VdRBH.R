setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(plyr)

Vd_only_genes <- read.table("../results/orthogroups/no_orthgroup_Vd_unmatched_VcGenes.txt", header = F)
annot <- read.csv("../results/annotations_GO.tsv", header = T,sep = "\t")
head(annot)
annot$Query.Sequence <- gsub("ID=","",annot$Query.Sequence)

Vd_only_annot <- merge(Vd_only_genes, annot, by.x="V1",by.y="Query.Sequence", all.x=T)

Vd2Vc <- read.table("../results/orthogroups/Vdarrrowii2Vc.blast", header = F, sep = "\t")
Vc2Vd <- read.table("../results/orthogroups/Vcorymbosum2Vd.blast", header = F, sep = "\t")
#
length(unique(Vd2Vc$V1))
# 32514 Vd genes matched to Vc 
length(unique(Vc2Vd$V1))
# 29493 Vc genes matched to Vd


Vd_spec <- read.table("../results/orthogroups/Vd_spec_genes.txt", header = F)
VdnotVc <- Vd_spec[!(Vd_spec$V1 %in% Vd2Vc$V1),]
VdnotVc <- VdnotVc[grep("ID=g",VdnotVc)]
VdnotVc_filter <- rea

Vd_only_annot_notVc <- Vd_only_annot[Vd_only_annot$V1 %in% gsub("ID=","",VdnotVc),]
