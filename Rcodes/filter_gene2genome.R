setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(plyr)

# total input sequences is 4502
blast <- read.table("../results/orthogroups/VdGene2Genome.blast", header = F, sep = "\t")
length(unique(blast$qseqid))
names(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
blast_filter <- blast[blast$length > 150, ]

order_evalue <- blast_filter[order(blast_filter$evalue, decreasing = F),]
order_evalue$pair <- paste0(order_evalue$qseqid,"_",order_evalue$sseqid)

uniq_blast <- unique(order_evalue[,c(1:2,13)])
uniq_blast <- join(uniq_blast, order_evalue[,3:13], by="pair", match="first")
length(unique(uniq_blast$qseqid))

# load RBH
rbh <- read.table("../results/reciprocal_best_hits.txt", header = F, sep = "\t")
match_Genes <- uniq_blast[uniq_blast$qseqid %in% rbh$V2,]
length(unique(match_Genes$qseqid))
unmatch_Vd_genes <- uniq_blast[!(uniq_blast$qseqid %in% rbh$V2),]
length(unique(unmatch_Vd_genes$qseqid))