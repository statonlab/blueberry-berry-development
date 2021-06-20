setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

filenames <- list.files("../results", pattern = "*node*.csv", full.names = T)
ldf <- lapply(filenames, function(x) read.csv(x, header = T))
# fix colnames, remove the sample level from column names
for (i in 1:7) {
  colnames(ldf[[i]]) <- gsub("_.*","",colnames(ldf[[i]]))
}
# combine all GO terms from the enrichment files
GOterms <- do.call(rbind, ldf)[,1:3]  
GOterms <- GOterms[!duplicated(GOterms$description),] #remove duplicated terms
# select the top five enriched GO terms from each sample
MEs <- gsub("../results/","",filenames)
MEs <- gsub(" default node.csv","",MEs)
BP <- GOterms[,c(1,3)]
for (i in 1:7) {
  sig <- ldf[[i]][order(ldf[[i]]$adjustedPValue, decreasing = F),]
  BP <- merge(BP, sig[1:5,c(2,3)], by="description", all.x=TRUE)
  names(BP)[i+2] <- MEs[i]
}

BP[is.na(BP)] <- 1
filterBP <- BP[rowSums(BP[,-c(1:2)]) < 7, ]
names(filterBP)[3:9] <- c("EPvLP_down","GvEP_up","GvLP_down","GvLP_up","GvR_down","GvR_up","LPvR_down")
filterBP <- filterBP[c("description","SUID","GvEP_up","GvLP_up","GvR_up","GvLP_down","GvR_down","EPvLP_down","LPvR_down")]
filterBP <- filterBP[order(filterBP$GvLP_up, decreasing = F),]
filterBP$description <- factor(filterBP$description)

# count how many GO terms in each file
BP_Total <- as.data.frame(matrix(NA,nrow=16,ncol=1))
rownames(BP_Total) <- MEs
for (i in 1:16) {
  c <- length(which(BP[,i+3]<0.05))
  BP_Total[i,1] <-c
}

# plot heatmap for top GO terms, using adjusted pvalues as colorscale
melt_BP_filtered <- melt(data = filterBP, id.vars = "description", measure.vars = names(filterBP)[3:9])
melt_BP_filtered$variable <- factor(melt_BP_filtered$variable,levels = c("GvEP_up","GvLP_up","GvR_up","GvLP_down","GvR_down","EPvLP_down","LPvR_down"))
melt_BP_filtered$description <- factor(melt_BP_filtered$description, levels = rev(filterBP$description))
p1 <- ggplot(melt_BP_filtered,aes(x=variable,y=description,fill=value))+
  geom_tile()+
  labs(fill='FDR') +
  scale_fill_gradient(low = "red", high = "white", limits = c(0,1)) +
  labs(x = "",y = "GO Term") +
  theme(axis.text.x = element_text(angle=75,hjust=1,vjust=1.0, size = 16),
        axis.text.y = element_text(size = 14))
ggsave("../img/topGO_all.pdf", width = 10, height = 5.6)
