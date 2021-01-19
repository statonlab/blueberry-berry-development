if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
############## install only needs to do once #########################
library(pathview)

# sample data
data(gse16873.d)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")
data(demo.paths)
data(paths.hsa)
head(paths.hsa,3)
i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                  pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = T, match.data = F, multi.state = T, same.layer = T)

##### my data ############################
# DEGs files
GvE <- read.csv("../results/DE_GvsEP.all.csv", header = T)
GvL <- read.csv("../results/DE_GvsLP.all.csv", header = T)
GvR <- read.csv("../results/DE_GvsR.all.csv", header = T)

# combine log2 fold change
DEG_data <- cbind(GvE[,c(1,3)], GvL[,3],GvR[,3])
# read the KEGG ortholog annotation
K_hits <- read.csv("../results/Vd1.2_Kgeneid.txt",na.strings = "", header = F, sep = "\t")
K_hits$V1 <- gsub("ID=","",K_hits$V1)

Kgenes <- merge(DEG_data, K_hits[,1:2], by.x = "X",by.y="V1", all.x=T)

K_data <- na.omit(Kgenes[,-c(1)])
K_uniq <- K_data[!duplicated(K_data$V2),]
rownames(K_uniq) <- K_uniq$V2
K_uniq <- as.matrix(K_uniq[,1:3])


pv.out <- pathview(gene.data = K_uniq, pathway.id = "00941",
                   species = "ko", out.suffix = "ko941")

pv.out <- pathview(gene.data = K_uniq, pathway.id = "00942",
                   species = "ko", out.suffix = "ko942")

pv.out <- pathview(gene.data = K_uniq, pathway.id = "00195",
                   species = "ko", out.suffix = "ko195")

pv.out <- pathview(gene.data = K_uniq, pathway.id = "00196",
                   species = "ko", out.suffix = "ko196")
