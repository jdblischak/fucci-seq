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

# select external validation samples
set.seed(99)
nvalid <- round(ncol(log2cpm.quant)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant), ii.valid)

pdata.nonvalid <- pdata[ii.nonvalid,]
pdata.valid <- pdata[ii.valid,]

#all.equal(rownames(pdata.valid), colnames(log2cpm.quant.valid))
log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
log2cpm.quant.valid <- log2cpm.quant[,ii.valid]
theta <- pdata$theta
names(theta) <- rownames(pdata)

log2cpm.nonvalid <- log2cpm.all[,ii.nonvalid]
log2cpm.valid <- log2cpm.all[,ii.valid]

# theta.nonvalid <- theta_moved[ii.nonvalid]
theta.nonvalid <- theta[ii.nonvalid]
theta.valid <- theta[ii.valid]

sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.101.rds")

log2cpm.quant.valid[rownames(log2cpm.quant.valid) %in% rownames(sig.genes), ]

### Scenario: fucci-theta, top 101 genes
data_training <- list(theta.nonvalid=theta.nonvalid,
                      log2cpm.quant.nonvalid=log2cpm.quant.nonvalid,
                      log2cpm.nonvalid=log2cpm.nonvalid,
                      pdata.nonvalid=pdata.nonvalid,
                      fdata=fdata)

data_withheld <- list(theta.valid=theta.valid,
                      log2cpm.quant.valid=log2cpm.quant.valid,
                      log2cpm.valid=log2cpm.valid,
                      pdata.valid=pdata.valid,
                      fdata=fdata)

saveRDS(data_training, file="data/results/data_training.rds")
saveRDS(data_withheld, file="data/results/data_withheld.rds")





############# <- get training partitions
# get predicted times
# set training samples
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")
folds <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
fold_indices <- folds$partitions

saveRDS(fold_indices, file="data/results/fold_indices.rds")
