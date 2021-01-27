setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(ggplot2)

## load data
gene_pair <- read.csv("../results/DupGenes.csv", header = T)
names(gene_pair)
ggplot(data = gene_pair, aes(x=Type, y=Gene.pairs, fill=Type))+
  geom_bar(stat = "identity")+
  labs(title = "", y = "Number of gene pairs", x = "Type")+
  geom_text(aes(y =Gene.pairs+1000, label = Gene.pairs), color = "black", size = 5)+
  theme_classic(base_size = 16)+
  theme(legend.position="none")
ggsave("../img/GeneDup.pdf")
