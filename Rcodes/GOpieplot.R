setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(ggplot2)
library(RColorBrewer)

GO_table <- read.csv("../results/annotations_GO.tsv", header = T, sep = "\t")
s <- strsplit(GO_table$EggNOG.GO.Biological, split = ",") # split the GO terms
anno.BP <- data.frame(locus = rep(GO_table$Query.Sequence, sapply(s, length)), GO = unlist(s))
L1BP <- anno.BP[grep("(L=1)", anno.BP$GO, fixed=TRUE),]

s <- strsplit(GO_table$EggNOG.GO.Cellular, split = ",") # split the GO terms
anno.CC <- data.frame(locus = rep(GO_table$Query.Sequence, sapply(s, length)), GO = unlist(s))
L1CC <- anno.CC[grep("(L=1)", anno.CC$GO, fixed=TRUE),]

s <- strsplit(GO_table$EggNOG.GO.Molecular, split = ",") # split the GO terms
anno.MF <- data.frame(locus = rep(GO_table$Query.Sequence, sapply(s, length)), GO = unlist(s))
L1MF <- anno.MF[grep("(L=1)", anno.MF$GO, fixed=TRUE),]

# plot pie chart
pieplot <- function(data, plotTitle,outpath){
df.data <- data.frame(table(data$GO))
df.data$Freq <- as.numeric(df.data$Freq)
df.data$pct <- round(df.data$Freq/sum(df.data$Freq), 4) * 100
df.data$Var1 <- gsub("GO:\\d+\\-|\\(L=1)","",df.data$Var1)
df.data <- df.data[order(df.data$Freq),]
df.data$Var1 <- paste0(df.data$Var1," (",df.data$pct,"%)")
df.data$Var1 <- factor(df.data$Var1, levels = df.data$Var1)

indplot<-ggplot(df.data,aes(x = "",y=Freq, fill = Var1))+ geom_bar(width = 1,stat="identity")+  
  coord_polar(theta = "y")+  
  labs(title= plotTitle, x="", y="")+
  theme_void()+ 
  scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(df.data$Var1)))+
  theme(legend.position="right",legend.title=element_blank(),
        plot.title = element_text(size = 16,hjust = 0.5),
        strip.text = element_text(size = 10), legend.text=element_text(size=10))
indplot
ggsave(outpath, width = 8, height = 5)
}

pieplot(L1BP, "Biological processes","../img/BPpiechart.pdf")  
pieplot(L1CC, "Cellular components","../img/CCpiechart.pdf")  
pieplot(L1MF, "Molecular functions","../img/MFpiechart.pdf")  

