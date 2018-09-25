#' @title Seurat code for cell cycle scoring
#' @description Compute phase specific cell cycle scores as outlined in Seurat

#' @param Y log2 normalized gene expression matrix
#'
#' @export
run_oscope <- function(Y, s.genes, g2m.genes, n.bin=25,
                       seed.use=1, random.seed=1) {

  library(Oscope)
  data_withheld <-readRDS("data/results/data_withheld.rds")
  df <- readRDS("data/eset-filtered.rds")
  counts <- exprs(df)
  samples_in_withheld <- match(colnames(data_withheld$log2cpm.quant.valid),colnames(counts))
  genes_in_withheld <- match(rownames(data_withheld$log2cpm.quant.valid), rownames(counts))
  counts_withheld <- counts[genes_in_withheld,samples_in_withheld]

  Sizes <- MedianNorm(counts_withheld)
  DataNorm <- GetNormalizedMat(counts_withheld, Sizes)
  # consdier CPM, mark genes for which the CPM is greater than 1
  # and a subset of these for the the observed log10 variance is greater than
  # the expected variation given the log10 gene expression mean
  MV <- CalcMV(Data = counts_withheld, Sizes = Sizes, MeanCutLow = 1)

  DataSubset <- DataNorm[MV$GeneToUse,]

  # check of the subsetted genes are all in the oscope genes
  #oscope_genes <- readRDS("data/cellcycle-genes-previous-studies/rds/leng-2015.rds")
  #sum(rownames(DataSubset) %in% oscope_genes$ensembl)
  #sum(rownames(counts_withheld) %in% oscope_genes$ensembl)

  DataInput <- NormForSine(DataSubset)
  SineRes <- OscopeSine(DataInput, parallel=TRUE)

  KMRes <- OscopeKM(SineRes, maxK = 10)

  ToRM <- FlagCluster(SineRes, KMRes, DataInput)
  print(ToRM$FlagID_bysine)
  print(ToRM$FlagID_byphase)
  print(ToRM$FlagID)
  KMResUse <- KMRes[-ToRM$FlagID]

  ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput, NCThre = 100)

  # obtain reordered dataset
  DataNorm2 <- DataNorm[,ENIRes[["cluster1"]]]

  fdata <- fData(df)
  ensembl <- rownames(fdata)
  genes_selected <- c("CDK1", "TOP2A", "UBE2C", "HIST1H4C")

  genes_ensembl <- ensembl[fdata$name %in% genes_selected]
  par(mfrow=c(2,2))
  for (g in 1:length(genes_ensembl)) {
    plot(DataNorm2[which(rownames(DataNorm2) == genes_ensembl[g]),])
  }

}



# use oscope to recover cell times
# then fit trendfilter, then order them by cyclical patterns...


