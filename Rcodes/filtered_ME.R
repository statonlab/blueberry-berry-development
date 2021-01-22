setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

# DEGs files
GvE <- read.csv("../results/DE_GvsEP.csv", header = T)
GvL <- read.csv("../results/DE_GvsLP.csv", header = T)
GvR <- read.csv("../results/DE_GvsR.csv", header = T)
EvL <- read.csv("../results/DE_EvsL.csv", header = T)
EvR <- read.csv("../results/DE_EPvsR.csv", header = T)
LvR <- read.csv("../results/DE_LvsR.csv", header = T)

all_DE <- c(as.character(GvE$X), as.character(GvL$X), as.character(GvR$X), as.character(EvL$X), as.character(EvR$X), as.character(LvR$X))

# load network
ME11 <- read.table("../results/CytoscapeInput-edges-ME11.txt", header = T, sep = "\t")
filtered_ME11 <- ME11[ME11$fromNode %in% all_DE,]
write.table(filtered_ME11, "../results/Cytoscape_ME11_filtered.txt", quote = F, row.names = F, sep = "\t")

ME7 <- read.table("../results/CytoscapeInput-edges-ME7.txt", header = T, sep = "\t")
filtered_ME7 <- ME7[ME7$fromNode %in% all_DE,]
write.table(filtered_ME7, "../results/Cytoscape_ME7_filtered.txt", quote = F, row.names = F, sep = "\t")

ME8 <- read.table("../results/CytoscapeInput-edges-ME8.txt", header = T, sep = "\t")
filtered_ME8 <- ME8[ME8$fromNode %in% all_DE,]
write.table(filtered_ME8, "../results/Cytoscape_ME8_filtered.txt", quote = F, row.names = F, sep = "\t")

ME9 <- read.table("../results/CytoscapeInput-edges-ME9.txt", header = T, sep = "\t")
filtered_ME9 <- ME9[ME9$fromNode %in% all_DE,]
write.table(filtered_ME9, "../results/Cytoscape_ME9_filtered.txt", quote = F, row.names = F, sep = "\t")

ME14 <- read.table("../results/CytoscapeInput-edges-ME14.txt", header = T, sep = "\t")
filtered_ME14 <- ME14[ME14$fromNode %in% all_DE,]
write.table(filtered_ME14, "../results/Cytoscape_ME14_filtered.txt", quote = F, row.names = F, sep = "\t")

length(unique(filtered_ME14$fromNode))
length(unique(ME11$fromNode))
