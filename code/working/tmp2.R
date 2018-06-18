df <- readRDS(file="../data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

macosko <- readRDS("../data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")

log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("../output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
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
theta.nonvalid <- theta[ii.nonvalid]
theta.valid <- theta[ii.valid]


sig.genes <- readRDS("../output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.rds")
expr.sig <- log2cpm.quant.nonvalid[rownames(log2cpm.quant.nonvalid) %in% rownames(sig.genes)[1:10], ]

# set training samples
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")
parts <- partitionSamples(1:ncol(log2cpm.quant.nonvalid), runs=5,
                          nsize.each = rep(151,5))
part_indices <- parts$partitions




######## get training samples
source("peco/R/fit.trendfilter.generic.R")
source("peco/R/cycle.npreg.R")
source("code/utility.R")

Y_train <- expr.sig[,part_indices[[1]]$train]
theta_train <- theta.nonvalid[part_indices[[1]]$train]

# one iteration to update training sample labels
fit.train.nobin <- cycle.npreg.insample(Y = Y_train,
                                        theta = theta_train,
                                        nbins = NULL, ncores=10)

plot(theta_train,
     fit.train.nobin$cell_times_est[match(names(theta_train),
                                          names(fit.train.nobin$cell_times_est))])


# per-cell misclassification
labels_to_test <- names(fit.train.nobin$cell_times_est)
labels_ref <- names(theta_train)
dist_mat_ref <- circ.dist.neighbors(labels_ref, k=10)
dist_mat_to_test <- circ.dist.neighbors(labels_to_test, k=10)
ii.match <- match(colnames(dist_mat_ref), colnames(dist_mat_to_test))
dist_mat_to_test <- dist_mat_to_test[ii.match, ii.match]

N <- length(labels_ref)
dist_enrich <- sapply(1:N, function(i) {
  lab_self <- rownames(dist_mat_ref)[i]
  labs_ref <- rownames(dist_mat_ref)[which(dist_mat_ref[i,]==1)]
  labs_to_test <- rownames(dist_mat_to_test)[which(dist_mat_to_test[i,]==1)]
  #      labs_to_test_nonneighbors <- setdiff(rownames(dist_mat_to_test), c(labs_to_test, lab_self))
  #prop_neighbors_in_ref <- sum(labs_to_test %in% labs_ref)/length(labs_to_test)
  #prop_nonneighbors_in_ref <- sum(labs_to_test_nonneighbors %in% labs_ref)/length(labs_to_test_nonneighbors)
  #      return(prop_neighbors_in_ref/prop_nonneighbors_in_ref)
  TP <- sum(dist_mat_ref[i,]==1 & dist_mat_to_test[i,]==1)
  FP <- sum(dist_mat_ref[i,]==0 & dist_mat_to_test[i,]==1)
  FN <- sum(dist_mat_ref[i,]==1 & dist_mat_to_test[i,]==0)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  if ((precision+recall)==0) {
    F1score <- 0
  } else {
    F1score <- 2*precision*recall/(precision+recall)
  }
  return(F1score)
})
names(dist_enrich) <- rownames(dist_mat_ref)
summary(dist_enrich)



str(theta_train)




#############

Y_test <- expr.sig[,part_indices[[1]]$test]
theta_test <- theta.nonvalid[part_indices[[1]]$test]

fit.test.bytrain.fucci.nobin <- cycle.npreg.outsample(Y_test,
                                                      theta_est=fit.train.nobin$cell_times_est,
                                                      mu_est=fit.train.nobin$mu_est,
                                                      sigma_est=fit.train.nobin$sigma_est)

plot(theta_test,
     fit.test.bytrain.fucci.nobin$cell_times_est[match(names(theta_test),
                                                       names(fit.test.bytrain.fucci.nobin$cell_times_est))])

# per-cell misclassification
labels_to_test <- names(fit.test.bytrain.fucci.nobin$cell_times_est)
labels_ref <- names(theta_test)
dist_mat_ref <- circ.dist.neighbors(labels_ref, k=10)
dist_mat_to_test <- circ.dist.neighbors(labels_to_test, k=10)
ii.match <- match(colnames(dist_mat_ref), colnames(dist_mat_to_test))
dist_mat_to_test <- dist_mat_to_test[ii.match, ii.match]

N <- length(labels_ref)
dist_enrich <- sapply(1:N, function(i) {
  lab_self <- rownames(dist_mat_ref)[i]
  labs_ref <- rownames(dist_mat_ref)[which(dist_mat_ref[i,]==1)]
  labs_to_test <- rownames(dist_mat_to_test)[which(dist_mat_to_test[i,]==1)]
  #      labs_to_test_nonneighbors <- setdiff(rownames(dist_mat_to_test), c(labs_to_test, lab_self))
  #prop_neighbors_in_ref <- sum(labs_to_test %in% labs_ref)/length(labs_to_test)
  #prop_nonneighbors_in_ref <- sum(labs_to_test_nonneighbors %in% labs_ref)/length(labs_to_test_nonneighbors)
  #      return(prop_neighbors_in_ref/prop_nonneighbors_in_ref)
  TP <- sum(dist_mat_ref[i,]==1 & dist_mat_to_test[i,]==1)
  FP <- sum(dist_mat_ref[i,]==0 & dist_mat_to_test[i,]==1)
  FN <- sum(dist_mat_ref[i,]==1 & dist_mat_to_test[i,]==0)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  if ((precision+recall)==0) {
    F1score <- 0
  } else {
    F1score <- 2*precision*recall/(precision+recall)
  }
  return(F1score)
})
names(dist_enrich) <- rownames(dist_mat_ref)
summary(dist_enrich)


