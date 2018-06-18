# try out prediction using ranks
# problem: there's higher error toward the end;
# this is probabiliy due to our smooth fitting that results in similar beginngin with the end
# solution
# try out doing different cutoffs for rotations to optimize the likelihood

library(Biobase)
df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

#macosko <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")

log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
# import previously identifid cell cycle genes
# cyclegenes <- readRDS("../output/npreg-methods.Rmd/cyclegenes.rds")
# cyclegenes.names <- colnames(cyclegenes)[2:6]

# select external validation samples
set.seed(99)
nvalid <- round(ncol(log2cpm.quant)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant), ii.valid)

log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
log2cpm.quant.valid <- log2cpm.quant[,ii.valid]
theta <- pdata$theta
names(theta) <- rownames(pdata)

# theta.nonvalid <- theta_moved[ii.nonvalid]
theta.nonvalid <- theta[ii.nonvalid]
theta.valid <- theta[ii.valid]

sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.rds")
expr.sig <- log2cpm.quant.nonvalid[rownames(log2cpm.quant.nonvalid) %in% rownames(sig.genes)[1:10], ]



# get predicted times
# set training samples
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")
parts <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
part_indices <- parts$partitions


#source("peco/R/fit.trendfilter.generic.R")
#source("peco/R/cycle.npreg.update.R")
source("peco/R/cycle.npcirc.R")
source("code/utility.R")

# fitting training data
Y_train <- expr.sig[,part_indices[[2]]$train]
theta_train <- theta.nonvalid[part_indices[[2]]$train]
fit.train.nobin <- cycle.npreg.insample(Y = Y_train,
                                        theta = theta_train,
                                        ncores=10,
                                        insample=T)

# there are three obvious three chunks
# also the beginning meet the end
par(mfrow=c(2,5))
for (i in 1:10) {
  plot(x=fit.train.nobin$cell_times_est, fit.train.nobin$Y_ordered[i,])
  points(x=fit.train.nobin$cell_times_est, fit.train.nobin$mu_est[i,], col = "blue", cex=.6, pch=16)
}



# fitting test data
Y_test <- expr.sig[,part_indices[[2]]$test]
theta_test <- theta.nonvalid[part_indices[[2]]$test]
fit.test.bytrain.nobin <- cycle.npreg.outsample(Y_test=Y_test,
                                                sigma_est=fit.train.nobin$sigma_est,
                                                funs_est=fit.train.nobin$funs_est,
                                                ncores=12,
                                                maxiter=20,
                                                tol=1, verbose=TRUE)

par(mfrow=c(2,5))
for (i in 1:10) {
  plot(x=fit.test.bytrain.nobin$cell_times_est,
       y=fit.test.bytrain.nobin$Y_ordered[i,])
  points(x=fit.test.bytrain.nobin$cell_times_est,
         y=fit.test.bytrain.nobin$mu_est[i,], col = "blue", cex=.6, pch=16)
}




