dir <- "/project2/gilad/joycehsiao/fucci-seq/data/rnaseq-previous-studies/bottcher"

library(data.table)
df <- fread(file.path(dir,"Fucci_GeneExpression_gencode22_TPM.txt"))
header <- readLines(file.path(dir,"Fucci_GeneExpression_gencode22_TPM.txt"), n=1)
cell_id <- strsplit(header, split="\t", fixed=T)[[1]]

# remove transcript indices
genenames <- sapply(1:nrow(df), function(g) strsplit(df$V1[g], split=".", fixed=TRUE)[[1]][1])
counts <- df[,-1]
colnames(counts) <- cell_id

pdata <- fread(file.path(dir, "FUCCI_cell_experimental metadata.txt"))
pdata_filtered <- pdata[pdata$Discard==F,]
counts_filtered <- counts[which(rownames(counts) %in% pdata_filtered$cell_id)]

log2cpm_all <- t(log2(1+(10^6)*(t(counts)/colSums(counts))))
log2cpm_filtered <- log2cpm_all[,which(colnames(log2cpm_all) %in% pdata_filtered$cell_id)]

gene_detection <- rowSums(log2cpm_filtered==0)/ncol(log2cpm_filtered)
summary(gene_detection)
genes_to_include <- which(gene_detection > .1)
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


saveRDS(log2cpm_filtered, file.path(dir, "log2cpm_filtered.rds"))
saveRDS(log2cpm_quant, file.path(dir, "log2cpm_quant.rds"))
saveRDS(pdata_filtered, file.path(dir, "pdata_filtered.rds"))






