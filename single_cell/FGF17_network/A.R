library(Seurat)
name=readRDS("FGFprocessed.RDS")
fgf17 <- subset(name, sample%in%"FGF17")
write.table(t(as.matrix(fgf17@assays$RNA@counts)), 
            paste0('network_FGF17/A_preparation/output/counts.tsv'), 
            sep = '\t', col.names = NA)

sink("network_FGF17/A_preparation/output/sessionInfo.txt")
sessionInfo()
sink()