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

theta.nonvalid <- theta[ii.nonvalid]
theta.valid <- theta[ii.valid]

sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.476.rds")

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



############ <- Make samples without noisy labels
library(Biobase)
df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

counts <- counts[,order(pdata$theta)]
log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")

# first check the theta in pdata
pca <- prcomp(cbind(pdata$gfp.median.log10sum.adjust,
                    pdata$rfp.median.log10sum.adjust), scale=TRUE)
pca_df <- cbind(pca$x[,1],pca$x[,2])
rownames(pca_df) <- rownames(pdata)

library(circular)
theta_check <- as.numeric(coord2rad(pca_df))
theta_check <- 2*pi-theta_check
#plot(theta_check, pdata$theta)
names(theta_check) <- rownames(pdata)

dist_to_origin <- sqrt(pca_df[,1]^2+pca_df[,2]^2)
which_out <- rownames(pdata)[which(scale(dist_to_origin) < -1)]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")

log2cpm.quant.filt <- log2cpm.quant[,!(colnames(log2cpm.quant) %in% which_out)]
#intersect(colnames(log2cpm.quant.filt), which_out)
log2cpm.all.filt <- log2cpm.all[,!(colnames(log2cpm.all) %in% which_out)]
#intersect(colnames(log2cpm.all.filt), which_out)
#all.equal(colnames(log2cpm.all.filt), colnames(log2cpm.quant.filt))
pdata.filt <- pdata[!(rownames(pdata) %in% which_out),]
#all.equal(rownames(pdata.filt), colnames(log2cpm.quant.filt))

theta <- theta_check
names(theta) <- rownames(pdata)
theta.filt <- theta[!(names(theta) %in% which_out) ]
#all.equal(names(theta.filt), rownames(pdata.filt))


# select external validation samples
set.seed(99)
nvalid <- round(ncol(log2cpm.quant.filt)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant.filt), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant.filt), ii.valid)

pdata.filt.nonvalid <- pdata.filt[ii.nonvalid,]
pdata.filt.valid <- pdata.filt[ii.valid,]

theta.filt.nonvalid <- theta.filt[ii.nonvalid]
theta.filt.valid <- theta.filt[ii.valid]

#all.equal(rownames(pdata.valid), colnames(log2cpm.quant.valid))
log2cpm.quant.filt.nonvalid <- log2cpm.quant.filt[,ii.nonvalid]
log2cpm.quant.filt.valid <- log2cpm.quant.filt[,ii.valid]

log2cpm.filt.nonvalid <- log2cpm.all.filt[,ii.nonvalid]
log2cpm.filt.valid <- log2cpm.all.filt[,ii.valid]


sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.476.rds")

data_filt_training <- list(theta.nonvalid=theta.filt.nonvalid,
                      log2cpm.quant.nonvalid=log2cpm.quant.filt.nonvalid,
                      log2cpm.nonvalid=log2cpm.filt.nonvalid,
                      pdata.nonvalid=pdata.filt.nonvalid,
                      fdata=fdata)

data_filt_withheld <- list(theta.valid=theta.filt.valid,
                      log2cpm.quant.valid=log2cpm.quant.filt.valid,
                      log2cpm.valid=log2cpm.filt.valid,
                      pdata.valid=pdata.filt.valid,
                      fdata=fdata)

saveRDS(data_filt_training, file="data/results/data_filt_training.rds")
saveRDS(data_filt_withheld, file="data/results/data_filt_withheld.rds")





############# <- get training partitions
# get predicted times
# set training samples
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")
folds <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
fold_indices <- folds$partitions

saveRDS(fold_indices, file="data/results/fold_indices.rds")
