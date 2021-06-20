setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(DESeq2)
library(ggplot2)
library(Rmisc)

# Load HTseq count data
mydata <- read.table("../counts/gene_counts_Vdarrowii.txt", header = TRUE, row.names = 1)
head(mydata)

# Import data description table
MydataDesign <- read.csv("../dataInfo/sample_info.csv",check.names = F, header = T, stringsAsFactors = T, row.names = 1)
head(MydataDesign)
mydata <- mydata[,sort(colnames(mydata))]
colnames(mydata) <- gsub("_","",colnames(mydata))

MydataDesign <- MydataDesign[sort(row.names(MydataDesign)),]
row.names(MydataDesign)
rownames(MydataDesign) == colnames(mydata)
# remove space for samples
MydataDesign$Stage <- gsub(" ","_",MydataDesign$Stage)
MydataDesign$Condition <- gsub(" ","_",MydataDesign$Condition)

dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~Condition)
dim(dds)
dds <- dds[rowSums(counts(dds))>1, ]
dim(dds)


# PCA plot
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = "Stage", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Stage)) +
  geom_point(size=3) + #geom_text(aes(label=name),hjust=0, vjust=0)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_discrete(breaks=c("Early green","Coloring","Ripening"))+
  theme_classic(base_size = 12)
ggsave("../img/PCA_vd.png")

x = "g12145"
data <- plotCounts(dds, gene=x,intgroup=c("Stage"), returnData=TRUE)
datasum <- summarySE(data, measurevar = "count", groupvars = "Stage" )
plot <-
  ggplot(datasum, aes(x=Stage, y=count, group=1)) +
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=count-se,ymax=count+se), width=0.1)+
  #scale_x_discrete(limits=c("Early green","Coloring","Ripening"))+
  labs(title = x, y = "normalized count")
plot

### differential expression analysis
# recode the sample conditions based on PCA. The coloring stage was redefined as early pink and late pink
# rerun dds
ddsMF = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~ Stage + Condition)
dim(ddsMF)
ddsMF <- ddsMF[rowSums(counts(ddsMF))>1, ]
dim(ddsMF)
levels(ddsMF$Condition)
# design(ddsMF) <- formula(~ Condition + Stage)
ddsMF <- DESeq(ddsMF, test = "LRT", reduced = ~ Condition)
resMF <- results(ddsMF)
head(resMF)
summary(resMF)
resultsNames(ddsMF)
dim(resMF[which(resMF$padj < 0.05),])
resultsNames(ddsMF)
resColor <- results(ddsMF, test = "Wald",contrast=c("Condition","Green","Early_pink"), lfcThreshold = 1, alpha = 0.05)
summary(resColor)

resColor <- results(ddsMF, test = "Wald",contrast=c("Condition","Green","Late_pink"), lfcThreshold = 1, alpha = 0.05)
summary(resColor)

resRipe <- results(ddsMF, test = "Wald",contrast=c("Condition","Green","Ripening"))
summary(resRipe)

## single factor wald test
dds <- DESeq(dds)
resdds <- results(dds)
summary(resdds)
resGreen <- results(dds, contrast = c("Condition","Early_pink","Green"), alpha = 0.05, lfcThreshold = 1)
summary(resGreen)
DE_green <- as.data.frame(resGreen[which(resGreen$padj < 0.05),])
DE_green <- as.data.frame(resGreen)
write.csv(DE_green,"../results/DE_GvsEP.all0315.csv")
resColor <- results(dds, contrast = c("Condition","Ripening","Early_pink"), lfcThreshold = 1, alpha = 0.05)
summary(resColor)
DE_Color <- as.data.frame(resColor[which(resColor$padj < 0.05),])
write.csv(DE_Color,"../results/DE_EPvsR.csv")

resRipe <- results(dds, contrast = c("Condition","Ripening","Green"), lfcThreshold = 1, alpha = 0.05)
summary(resRipe)
DE_Ripe <- as.data.frame(resRipe[which(resRipe$padj < 0.05),])
write.csv(DE_Ripe,"../results/DE_GvsR.csv")

resDE <- results(dds, contrast = c("Condition","Early_pink","Ripening"), lfcThreshold = 1, alpha = 0.05)
summary(resDE)
