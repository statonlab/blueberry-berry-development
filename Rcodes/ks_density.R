setwd("~/Desktop/jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

kaks <- read.table("../../Vd_wgd/Rsimsii.wgd.kaks", header = T, sep = "\t", stringsAsFactors = F)
# Kernel Density Plot
d <- density(na.omit(as.numeric(kaks[-1,]$Ks))) # returns the density data
plot(d, xlim=c(0,2), ylim=c(0,1)) # plots the results
