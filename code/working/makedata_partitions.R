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
nvalid <- round(ncol(log2cpm.quant)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant), ii.valid)

log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
log2cpm.quant.valid <- log2cpm.quant[,ii.valid]

folds <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
fold_indices <- folds$partitions

saveRDS(fold_indices, file="data/results/fold_indices.rds")



#######################----- validation sample of random indivdiual
for (ind in unique(pdata$chip_id)) {
  set.seed(99)
  #  nvalid <- round(ncol(log2cpm.quant)*.15)
  ii.valid <- c(1:nrow(pdata))[which(pdata$chip_id == ind)]
  ii.nonvalid <- c(1:nrow(pdata))[which(pdata$chip_id != ind)]

  pdata.nonvalid <- pdata[ii.nonvalid,]
  pdata.valid <- pdata[ii.valid,]

  # log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
  # log2cpm.quant.valid <- log2cpm.quant[,ii.valid]
  # theta <- pdata$theta
  # names(theta) <- rownames(pdata)

  # log2cpm.nonvalid <- log2cpm.all[,ii.nonvalid]
  # log2cpm.valid <- log2cpm.all[,ii.valid]
  #
  # theta.nonvalid <- theta[ii.nonvalid]
  # theta.valid <- theta[ii.valid]

  # #sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.476.rds")
  # data_training <- list(theta.nonvalid=theta.nonvalid,
  #                       log2cpm.quant.nonvalid=log2cpm.quant.nonvalid,
  #                       log2cpm.nonvalid=log2cpm.nonvalid,
  #                       pdata.nonvalid=pdata.nonvalid,
  #                       fdata=fdata)
  #
  # data_withheld <- list(theta.valid=theta.valid,
  #                       log2cpm.quant.valid=log2cpm.quant.valid,
  #                       log2cpm.valid=log2cpm.valid,
  #                       pdata.valid=pdata.valid,
  #                       fdata=fdata)
  #
  # saveRDS(data_training, file=paste0("data/results/ind_",ind,"_data_training.rds"))
  # saveRDS(data_withheld, file=paste0("data/results/ind_",ind,"_data_withheld.rds"))

  ############# <- get training partitions
  # split by individaul
  # # get predicted times
  # # set training samples
  # source("peco/R/primes.R")
  # source("peco/R/partitionSamples.R")
  # folds <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
  #                           nsize.each = c(rep(round(ncol(log2cpm.quant.nonvalid)/5),4),
  #                                          ncol(log2cpm.quant.nonvalid)-sum(rep(round(ncol(log2cpm.quant.nonvalid)/5),4))))
  fold_indices <- lapply(1:length(unique(pdata.nonvalid$chip_id)), function(i) {
    ind_test <- unique(pdata.nonvalid$chip_id)[i]
    test <- which(pdata.nonvalid$chip_id==ind_test)
    train <- which(pdata.nonvalid$chip_id!=ind_test)
    return(list(test=test, train=train))
  })
#  folds$partitions

  saveRDS(fold_indices, file=paste0("data/results/ind_",ind,"_fold_indices.rds"))
}




