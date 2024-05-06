library(Seurat)
name=readRDS("FGFprocessed.RDS")
fgf8 <- subset(name, sample%in%"FGF8")
write.table(t(as.matrix(fgf8@assays$RNA@counts)), 
            paste0('network_FGF8/A_preparation/output/counts.tsv'), 
            sep = '\t', col.names = NA)

sink("network_FGF8/A_preparation/output/sessionInfo.txt")
sessionInfo()
sink()