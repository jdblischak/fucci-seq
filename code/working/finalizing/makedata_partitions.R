################################################
# Description:
#   Partition samples to training and validation
################################################

library(Biobase)
df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

#macosko <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")
counts <- counts[,order(pdata$theta)]
log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")


#######################----- validation sample of random cells
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")

# select external validation samples
set.seed(99)
nsamples <- ncol(log2cpm.quant)
nsizes <- c(nsamples-148,148)
folds <- 6

seeds <- primes(folds)
out <- lapply(1:folds, function(f) {
  set.seed(seeds[f])
  indices <- sample(1:nsamples, size=nsizes, replace = F)
  indices <- sort(indices)
  part_indices <- list(train=indices,
                       test=setdiff(1:nsamples,indices))
  return(part_indices)
  })

saveRDS(out, file="data/results/mixed_overall_fold_indices.rds")



