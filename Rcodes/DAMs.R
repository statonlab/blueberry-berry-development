setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

DAM_blast <- read.table("../results/Vd_DAMs.tsv", header = F, sep = "\t")

# filter genes with 50% identity
DAM_blast_filter <- DAM_blast[DAM_blast$V4 > 50,]
length(unique(DAM_blast_filter$V1)) # 17 genes left

# try filter by a python script with 50% CIP and 70% CALP
DAM_filter <- read.table("../results/Vd_DAM_filtered.tsv", header = F, sep = "\t")
