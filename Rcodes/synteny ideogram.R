setwd("~/Desktop/Jiali/UTK/blueberry/scaffolds/")
require(RIdeogram)

synteny_links <- read.table("maps/Vc2Vd_Feb2021.coord.aligncoords.filtered", header = F, sep = "\t")
names(synteny_links) <- c("Species_1","Start_1","End_1","Species_2","Start_2","End_2","fill")
synteny_links$Species_1 <- gsub("chr","",synteny_links$Species_1)
karyotpe_1 <- read.table("circos/karyotype.chr1-12.txt", header = F, sep="\t")
names(karyotpe_1)[3:6] <- c("Chr","id","Start","End")
karyotype_2 <- read.table("circos/karyotype.scaff.txt", header = F, sep = "\t")
names(karyotype_2)[3:6] <- c("Chr","id","Start","End")
karyotype_2$Chr <- paste0("chr",1:12)

karyotype_dual <- rbind(karyotpe_1[9:12,c(3,5:6)],karyotype_2[9:12,c(3,5:6)])
karyotype_dual$fill <- rep("969696",8)
karyotype_dual$species <- c(rep("Vcory",4), rep("Vdarr",4))
karyotype_dual$size <- rep(12,8)
karyotype_dual$color <- rep("252525",8)
karyotype_dual$Chr <- c(9:12,9:12)

synteny_dual <- data.frame(matrix(ncol = 7, nrow = 0))
names(synteny_dual) <- names(synteny_dual_comparison)
for (i in 9:12) {
  links <- synteny_links[which(synteny_links$Species_1 == karyotpe_1$Chr[i] & synteny_links$Species_2 == karyotype_2$Chr[i]),]
  links$Species_1 <- i-8
  links$Species_2 <- i-8
  links$fill <- "cccccc"
  synteny_dual <- rbind(synteny_dual, links)
}

ideogram(karyotype = karyotype_dual, synteny = synteny_dual, width = 200)
convertSVG("chromosome.svg", device = "png")
