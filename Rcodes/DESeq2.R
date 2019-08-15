setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(DESeq2)
library(ggplot2)

# Load HTseq count data
mydata <- read.table("../counts/gene_counts_Vcorymbosum.txt", header = TRUE, row.names = 1)
head(mydata)

# Import data description table
MydataDesign <- read.csv("../dataInfo/sample_info.csv",check.names = F, header = T, stringsAsFactors = T, row.names = 1)
head(MydataDesign)
mydata <- mydata[,sort(colnames(mydata))]
colnames(mydata) <- gsub("_","",colnames(mydata))

MydataDesign <- MydataDesign[sort(row.names(MydataDesign)),,drop=F]
row.names(MydataDesign)
rownames(MydataDesign) == colnames(mydata)

dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~Stage)
dim(dds)
dds <- dds[rowSums(counts(dds))>1, ]
dim(dds)
dds <- DESeq(dds)

# PCA plot
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = "Stage", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Stage)) +
  geom_point(size=3) + geom_text(aes(label=name),hjust=0, vjust=0)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_discrete(breaks=c("Early green","Coloring","Ripening"))+
  theme_classic(base_size = 12)
