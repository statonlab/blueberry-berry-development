setwd("~/Desktop/Jiali/UTK/blueberry/blueberry-berry-development/Rcodes/")

target_table <- read.table("../results/miR/psRNATargetJob-1623903871905904.txt", sep = "\t", skip = 1, header = T)

filter_table <- target_table[target_table$Expectation < 3,]

length(unique(filter_table$Target_Acc.))
length(unique(filter_table$miRNA_Acc.))

length(unique(target_table$Target_Acc.))

Vd_spec <- read.table("../results/orthogroups/Vd_spec_unmatched_Vcgenome.txt", header = F)

GO_anno$Query.Sequence <- gsub("ID=","",GO_anno$Query.Sequence)
Vd_function <- merge(Vd_spec, GO_anno, by.x="V1", by.y = "Query.Sequence", all.x = T)

write.csv(Vd_function, "../results/orthogroups/Vd_spec_unmatched_Vcgenome_anno.csv", row.names = F)

write.csv(filter_table, "../results/miR/miRNA_targets.filtered.csv", row.names = F)
