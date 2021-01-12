setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")
library(ggplot2)
library(reshape2)

# plot number of DEG
data <- data.frame(samples = c("GvEP","GvLP","GvR","EPvLP","EPvR","LPvR"), 
                   down = c(7,1607,2528,619,1369,264),
                   up = c(134,751,1067,14,318,149))
data_melt <- melt(data)
data_melt$samples <- factor(data_melt$samples, levels = c("GvEP","GvLP","GvR","EPvLP","EPvR","LPvR"))

ggplot(data=data_melt, aes(x=samples, y=value, group=variable,fill=variable)) +
  geom_col(position="dodge") +
  labs(title = "",y="Number of DEGs", x = "")+
  scale_fill_manual(values=c("#4575B495","#D7302795"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1,size = 16),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 16))
ggsave("../img/DEG number.png")
