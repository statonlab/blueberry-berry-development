setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

GO_table <- read.csv("../results/annotations_GO.tsv", header = T, sep = "\t")
GO_table$Query.Sequence <- gsub("ID=","",GO_table$Query.Sequence)
s <- strsplit(GO_table$EggNOG.GO.Biological, split = ",") # split the GO terms
anno.BP <- data.frame(locus = rep(GO_table$Query.Sequence, sapply(s, length)), GO = unlist(s))
anno.BP$GO <- gsub("GO:|\\-.*","",anno.BP$GO)
anno.BP$GO <- as.numeric(anno.BP$GO)
write.table(anno.BP, file="../results/go_BP_anno.txt", sep = " = ", quote = F , row.names = F, col.names = F)
length(intersect(anno.BP$locus, rownames(mydata)))
