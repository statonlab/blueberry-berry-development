setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(DESeq2)
library(ggplot2)

novel_miR <- read.csv("../results/miR/result_27_01_2021_t_10_08_03.csv", header = T, stringsAsFactors = F)
filtered_miR <- novel_miR[which(novel_miR$miRDeep2_score > 3 & novel_miR$randfold_p.value == 'yes' & novel_miR$mature_count > 10),]
filtered_miR$length <- nchar(filtered_miR$consensus)

# keep plant miRs
organism_tree <- read.table("../dataInfo/miRbase organisms.txt", header = T, sep = "\t")
plant_ID <- organism_tree[grep("Viridiplantae", organism_tree$tree),]
filtered_miR$species <- sapply(strsplit(filtered_miR$example, split="-"), "[",1)
filtered_miR_plant <- filtered_miR[filtered_miR$species %in% plant_ID$organism,]

# filter by size
potential_miR <- filtered_miR_plant[which(filtered_miR_plant$length > 20 & filtered_miR_plant$length < 25),]
length(unique(potential_miR$example))
unmatch_miR <- filtered_miR[which(filtered_miR$species == ""),]
potential_miR <- rbind(potential_miR, unmatch_miR)

# size distribution
barplot(table(potential_miR$length))

# expression analysis
S10 <- read.csv("../results/miR/miRNAs_expressed_all_samples_1612286023_S10.csv", header = T, sep = "\t")
S14 <- read.csv("../results/miR/miRNAs_expressed_all_samples_1612286325_S14.csv", header = T, sep = "\t")
S16 <- read.csv("../results/miR/miRNAs_expressed_all_samples_1612286545_S16.csv", header = T, sep = "\t")
S17 <- read.csv("../results/miR/miRNAs_expressed_all_samples_1612286750_S17.csv", header = T, sep = "\t")

miR_data <- cbind(S10[,1:2],S14[,1:2],S16[,1:2],S17[,1:2])
miR_data <- miR_data[1:8985,c(1,2,4,6,8)]
names(miR_data) <- c("miRNA","S10","S14","S16","S17")
potential_miR_Data <- miR_data[miR_data$miRNA %in% potential_miR$id,]
rownames(potential_miR_Data) <- potential_miR_Data$miRNA

mySampleTable <- data.frame(row.names = colnames(potential_miR_Data[,2:5]), Stage = c("coloring",rep("ripening",3)))

dds_m = DESeqDataSetFromMatrix(countData = potential_miR_Data[,2:5], colData = mySampleTable, design = ~Stage)
dds_m = DESeq(dds_m)
dds_m = dds_m[rowSums(counts(dds_m))>1, ]
dds_m= DESeq(dds_m)
dim(dds_m)
rld <- rlog(dds_m, blind = FALSE)

rl.heatmap <- assay(rld)
library(pheatmap)

pheatmap(rl.heatmap, filename = "../img/miR_heatmap.pdf",cluster_cols = F,annotation_col=mySampleTable[,1, drop=F], 
         border_color = NA,show_rownames = F,color = colorRampPalette(c("white","red"))(50),
         width = 4, height = 8)

# accumulative expression
expression <- as.data.frame(rl.heatmap)
expression_data <- data.frame("x" = c(1:length(expression$S10)))
for (a in 1:4) {
  y <- cumsum(expression[,a]) 
  expression_data[,a+1] <- expression[,a]/sum(expression[,a])
}
expression_data$x <- rownames(expression)
# curve_data <- data.frame("x" = rownames(expression))
# curve_data$D3 <- rowMeans(expression_data[,2:4])
# curve_data <- curve_data[order(curve_data$D3, decreasing = T),]
# curve_data$x <- factor(curve_data$x, levels = curve_data$x)
# curve_data$cum <- cumsum(curve_data$D3)
# ggplot(curve_data[1:15,], aes(x= x, weight = D3))+
#   geom_bar(width = 0.8, fill = "#00B0F6") +
#   geom_line(aes(x=x, y = cum, group=1), colour="#00B0F6") +
#   theme_classic() +
#   theme(axis.text.x=element_text(angle=60, hjust=1),
#         axis.text = element_text(colour = "black"))+xlab("")+
#   geom_hline(yintercept =0.9, linetype="dashed", color="black")+
#   ylab("expression fraction")
# 
# curve_data$D7 <- rowMeans(expression_data[,5:7])
# curve_data$T1 <- rowMeans(expression_data[,8:10])
# curve_data$T2 <- rowMeans(expression_data[,11:13])
# curve_data$T3 <- rowMeans(expression_data[,14:15])
# ## color code: T3: "#F8766D" D7:"#A3A500" T2:"#00BF7D" D3:"#00B0F6" T1:"#E76BF3"

sort_curve <- expression_data[order(expression_data$V2, decreasing = T),]
sort_curve_cum <- data.frame("x" = c(1:length(expression_data$x)))
sort_curve_cum$S10 <- cumsum(sort_curve$V2)
sort_curve_cum$S14 <- cumsum(sort_curve$V3[order(sort_curve$V3, decreasing = T)])
sort_curve_cum$S16 <- cumsum(sort_curve$V4[order(sort_curve$V4, decreasing = T)])
sort_curve_cum$S17 <- cumsum(sort_curve$V5[order(sort_curve$V5, decreasing = T)])

library(reshape2)
melt_curve <- melt(sort_curve_cum, id.vars = "x")
melt_curve$variable <- factor(melt_curve$variable, levels = c("S10","S14","S16","S17"))

ggplot(data = melt_curve, aes(x=x, y = value, group=variable)) +
  geom_line( aes(color=variable), size =1)+
  theme_classic(base_size = 14)+
  xlab("Number of microRNA")+
  ylab("Cumulative expression fraction")+
  geom_hline(yintercept =0.9, linetype="dashed", color="black")+
  scale_color_manual(values = c("#E76BF3","#00BF7D","#F8766D","#00B0F6"))+
  guides(color = guide_legend(title = "Sample"))

ggsave("../img/allmiR_curve.pdf", height = 4, width = 4)

# output miRNA
miR_anno <- merge(expression, potential_miR, by.x = "row.names", by.y = "id", all.x=T)
write.csv(miR_anno, "../results/miR_rld_annot.csv")
