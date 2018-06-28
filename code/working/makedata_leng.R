df <- readRDS(file="data/rnaseq-previous-studies/leng/HumanLengESC.rds")
pdata <- pData(df)
counts <- exprs(df)

pdata_filtered <- pdata[pdata$cell_state != "H1",]
pdata_filtered$cell_state <- drop.levels(pdata_filtered$cell_state)

counts_filtered <- counts[,colnames(counts) %in% pdata_filtered$sample_id]
dim(counts_filtered)

log2cpm_filtered <- log2(10^6*t(t(counts_filtered)/colSums(counts_filtered))+1)

counts_ensg_filtered <- counts_filtered[grep("ERCC", rownames(counts_filtered), invert=TRUE),]

genes_to_include <- which(rowMeans(counts_ensg_filtered > 0) > .9)
log2cpm_filtered <- log2cpm_filtered[genes_to_include,]


library(parallel)
log2cpm_quant <- do.call(rbind, mclapply(1:nrow(log2cpm_filtered), function(g) {
  yy <- log2cpm_filtered[g,]
  is.zero <- which(yy == 0)
  qq.map <- qqnorm(yy)
  yy.qq <- qq.map$x
  yy.qq[is.zero] <- sample(qq.map$x[is.zero])
  return(yy.qq)
}, mc.cores=8) )
rownames(log2cpm_quant) <- rownames(log2cpm_filtered)

saveRDS(log2cpm_quant,
        file = "data/rnaseq-previous-studies/leng/log2cpm_quant.rds")
saveRDS(log2cpm_filtered,
        file = "data/rnaseq-previous-studies/leng/log2cpm_filtered.rds")
saveRDS(pdata_filtered,
        file = "data/rnaseq-previous-studies/leng/pdata_filtered.rds")
