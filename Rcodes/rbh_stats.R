setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

RBH <- read.table("../results/reciprocal_best_hits.txt", header = F)

length(unique(RBH$V1))
length(unique(RBH$V2))


RBH_Vc <- read.table("../results/reciprocal_best_hitsVc2Vd.txt", header = F)
RBH_Vc_p <- RBH_Vc[-grep("h2",RBH_Vc$V1),]
length(unique(RBH_Vc_p$V1))
length(unique(RBH_Vc_p$V2))

RBH_Vc_h <- RBH_Vc[grep("h2",RBH_Vc$V1),]
length(unique(RBH_Vc_h$V1))
length(unique(RBH_Vc_h$V2))
