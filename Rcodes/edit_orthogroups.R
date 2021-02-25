setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

orthgroups <- read.table("../results/orthogroups/Orthogroups.tsv", sep = "\t", header = T)

Vd_spec <- orthgroups[which(orthgroups$V_corymbosum_Draper_v1.0.proteins.nameTruncated == ""),]
s <- strsplit(as.character(Vd_spec$Vdarrowii_all_proteins), split = ",") # split the GO terms
genes <- data.frame(ID = unlist(s))
genes$ID <- gsub(" ","",genes$ID)

nonorthogenes <- read.table("../results/orthogroups/Orthogroups_UnassignedGenes.tsv", sep = "\t", header = T)
Vd_genes <- nonorthogenes[which(nonorthogenes$V_corymbosum_Draper_v1.0.proteins.nameTruncated == ""),]
names(Vd_genes)[3] <- "ID"

genes_all <- rbind(genes, Vd_genes[,3, drop=F])
genes_all$ID <- gsub(" ","",genes_all$ID)

write.table(genes_all,"../results/orthogroups/Vd_spec_genes.txt", row.names = F, col.names = F, quote = F)

# RBH 18907 genes


