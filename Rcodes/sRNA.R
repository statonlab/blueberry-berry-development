setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(plyr)
Result_table <- read.table("../results/hap2_sRNA.txt", header = F)
names(Result_table)[1:2] <- c("chr","Name")
blast_cluster <- read.table("../results/hap2_ShortStack_all_clusters_blastn0120.txt", header = F)
names(blast_cluster)[1:2] <- c("Name","ref")

blast_Result_table <- join(Result_table, blast_cluster[,1:2], by = "Name", match = "first")
length_dist <- table(Result_table$V12)
barplot(length_dist, main="All ShortStack clusters",
        xlab="length")

cluster_24bp <- blast_Result_table[grep("24",blast_Result_table$V12),]
cluster_24bp$category <- gsub("_.*","",cluster_24bp$ref)
length_24bp <- as.data.frame(table(cluster_24bp$category))
de<-data.frame("other",sum(is.na(cluster_24bp$category)))
names(de)<-c("Var1","Freq")
length_24bp <- rbind(length_24bp, de)
pie(length_24bp$Freq,labels = length_24bp$Var1, main="24-bp clusters category")
length_24bp$prop <- length_24bp$Freq/sum(length_24bp$Freq)
ggplot(length_24bp, aes(x = "", y = prop, fill = Var1)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  labs(fill='Category') +
  scale_fill_brewer(palette="Dark2") +
  theme_void()

cluster_20bp <- blast_Result_table[grep("20",blast_Result_table$V12),]
cluster_20bp$category <- gsub("_.*","",cluster_20bp$ref)
length_20bp <- as.data.frame(table(cluster_20bp$category))
de<-data.frame("other",sum(is.na(cluster_20bp$category)))
names(de)<-c("Var1","Freq")
length_20bp <- rbind(length_20bp, de)

cluster_21bp <- blast_Result_table[grep("21",blast_Result_table$V12),]
cluster_21bp$category <- gsub("_.*","",cluster_21bp$ref)
length_21bp <- as.data.frame(table(cluster_21bp$category))
de<-data.frame("other",sum(is.na(cluster_21bp$category)))
names(de)<-c("Var1","Freq")
length_21bp <- rbind(length_21bp, de)

cluster_22bp <- blast_Result_table[grep("22",blast_Result_table$V12),]
cluster_22bp$category <- gsub("_.*","",cluster_22bp$ref)
length_22bp <- as.data.frame(table(cluster_22bp$category))
de<-data.frame("other",sum(is.na(cluster_22bp$category)))
names(de)<-c("Var1","Freq")
length_22bp <- rbind(length_22bp, de)

cluster_23bp <- blast_Result_table[grep("23",blast_Result_table$V12),]
cluster_23bp$category <- gsub("_.*","",cluster_23bp$ref)
length_23bp <- as.data.frame(table(cluster_23bp$category))
de<-data.frame("other",sum(is.na(cluster_23bp$category)))
names(de)<-c("Var1","Freq")
length_23bp <- rbind(length_23bp, de)

cluster_N <- blast_Result_table[grep("N",blast_Result_table$V12),]
cluster_N$category <- gsub("_.*","",cluster_N$ref)
length_N <- as.data.frame(table(cluster_N$category))
de<-data.frame("other",sum(is.na(cluster_N$category)))
names(de)<-c("Var1","Freq")
length_N <- rbind(length_N, de)

names(length_20bp)[2] <- "20bp"
names(length_21bp)[2] <- "21bp"
names(length_22bp)[2] <- "22bp"
names(length_23bp)[2] <- "23bp"
names(length_24bp)[2] <- "24bp"
names(length_N)[2] <- "N"

clusters_cate <- join(length_20bp, length_21bp, by = "Var1", type = "full")
clusters_cate <- join(clusters_cate, length_22bp, by = "Var1", type = "full")
clusters_cate <- join(clusters_cate, length_23bp, by = "Var1", type = "full")
clusters_cate <- join(clusters_cate, length_24bp[,1:2], by = "Var1", type = "full")
clusters_cate <- join(clusters_cate, length_N, by = "Var1", type = "full")

library(reshape2)
data <- melt(clusters_cate)
data$Var1 <- factor(data$Var1, levels = c("miRNA","rRNA","snoRNA","splicing","tRNA","other","ta-siRNA"))
ggplot(data, aes(fill=Var1, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity") + theme_classic()+
  labs(x= "size", y = "Number of clusters", fill='Category') +
  scale_fill_brewer(palette="Dark2")
ggsave("../img/Hap2_sRNA_dist.pdf")

## Remove the known tRNA, rRNA, snoRNA and splicing, and invalid clusters
blast_Result_table$category <- gsub("_.*","",blast_Result_table$ref)
filtered_cluster <- blast_Result_table[!(blast_Result_table$category %in% c("tRNA","rRNA","snoRNA","splicing")),]
filtered_cluster_N <- filtered_cluster[which(filtered_cluster$V12 != "N"),]

## gene region
intragenic_site <- read.table("../results/hap2_ShortStack_gene_overlap.gff", sep = "\t", stringsAsFactors = F)
intragenic_ID <- sapply(strsplit(intragenic_site$V9,";"),"[",1)
intragenic_ID <- gsub("ID=","",intragenic_ID)
repeat_site <- read.table("../results/hap2_ShortStack_repeats.gff", sep = "\t", stringsAsFactors = F)
repeat_site$size <- sapply(strsplit(repeat_site$V9,";"), "[",2)
repeat_ID <- sapply(strsplit(repeat_site$V9,";"), "[",1)
repeat_ID <- gsub("ID=","",repeat_ID)

## clusters loci
table(filtered_cluster_N$Name %in% intragenic_ID)
table(filtered_cluster_N$Name %in% repeat_ID)
regions <- data.frame("Location" = c("intergenic","repeats","gene"), "Number" = c(188,174285, 27919))
regions$prob <- round((regions$Number/sum(regions$Number) * 100),1)

ggplot(regions, aes(x = "", y = prob, fill = Location)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0)+
  labs(fill='Location') +
  scale_fill_brewer(palette="Dark2") +
  #geom_text(aes(y =cumsum(prob) - 0.2*prob , label = paste0(prob,"%")), color = "black", size = 5)+
  theme_void(base_size = 14)
ggsave("../img/hap2_sRNA_pie.pdf")

