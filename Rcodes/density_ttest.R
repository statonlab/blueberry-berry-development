setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
#install.packages("vioplot")
library(vioplot)

# small RNA density
sRNA <- read.table("../dataInfo/ncRNA_all.bg", header = F, sep = "\t")
# repeats density
repeats <- read.table("../dataInfo/repeats_all.bg", header=F, sep = "\t")
# CG density 
GC <- read.table("../dataInfo/CG_chr.all.bg", header=F, sep="\t")

repeat_rich <- repeats[repeats$V4 > 60,]
GC_rich <- GC[GC$V4 > 46,]

repeat_low <- repeats[repeats$V4 < 29,]
GC_low <- GC[GC$V4 < 35,]

# statistical test
t.test(sRNA[repeats$V4>60,]$V4, sRNA[repeats$V4<29,]$V4)

vioplot(sRNA[repeats$V4>60,]$V4, sRNA[repeats$V4<29,]$V4, main="small RNA density",
        names=c("repeat-rich", "repeat-low"))

t.test(sRNA[GC$V4>45,]$V4, sRNA[GC$V4<35,]$V4)
vioplot(sRNA[GC$V4>45,]$V4, sRNA[GC$V4<35,]$V4, main="small RNA density",
        names=c("GC-rich", "GC-low"))
