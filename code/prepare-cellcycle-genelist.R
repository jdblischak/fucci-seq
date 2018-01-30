###### ---- Get Macosko 2015 genes from the human experiment
# Macosko identified 544 cell cycle genes using human cells
# the output below gives 558 ensembl gene IDs, because
# in several cases, a HGNC symbol is mapped to multiple ensembl gene ID
macosko <- xlsx::read.xlsx("data/cellcycle-genes-previous-studies/cellcyclegenes-macosko-2015.xlsx",
                           sheetName = 4, header = TRUE)
tmp <- vector("list", length(macosko))

for (i in 1:length(phases)) {
  macosko[[i]] <- as.character(macosko[[i]])
  macosko[[i]] <- gsub(" ", "", macosko[[i]])
  tmp[[i]] <- macosko[[i]][!is.na(macosko[[i]])]
}

tmp <- data.frame(hgnc = c(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]], tmp[[5]]),
                  phase = c( rep(names(macosko)[1], length(tmp[[1]])),
                             rep(names(macosko)[2], length(tmp[[2]])),
                             rep(names(macosko)[3], length(tmp[[3]])),
                             rep(names(macosko)[4], length(tmp[[4]])),
                             rep(names(macosko)[5], length(tmp[[5]])) ),
                  stringsAsFactors = FALSE)

library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

biomart <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "hgnc_symbol",
                 values = tmp$hgnc,
                 mart = ensembl)

anno <- merge(macosko_anno, biomart, by.x="hgnc", by.y="hgnc_symbol", all=TRUE)
anno$hgnc <- as.character(anno$hgnc)
anno$phase <- as.character(anno$phase)

colnames(anno) <- c("hgnc", "phase", "ensembl")

# change cell cycle phase labeling to be consistent with the labeling in Whitfield et al. 2002

anno$phase[anno$phase=="G1.S"] <- "G1/S"
anno$phase[anno$phase=="G2.M"] <- "G2"
anno$phase[anno$phase=="M.G1"] <- "M/G1"

saveRDS(anno, file = "data/cellcycle-genes-previous-studies/rds/macosko-2017.rds")
