#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}



library(circular)
library(conicfit)
library(Biobase)
library(dplyr)
library(matrixStats)

library(NPCirc)
library(smashr)
library(parallel)
library(genlasso)

## Load data
dir <- "/project2/gilad/joycehsiao/fucci-seq"
df <- readRDS(file.path(dir,"data/eset-final.rds"))
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

# import corrected intensities
pdata.adj <- readRDS(file.path(dir,"output/images-normalize-anova.Rmd/pdata.adj.rds"))

#log2cpm.detected <- log2cpm.all[rowMeans(log2cpm.all>0)>.8,]

source(file.path(dir,"code/utility.R"))

## cell time
theta <- readRDS(file.path(dir, "output/images-time-eval.Rmd/theta.rds"))

# ## order cells based on estimated time
partitions <- partitionSamples.n(c(1:length(theta)), runs=1,
                                 c(512,length(theta)-512))[[1]][[1]]

ord <- order(theta)
log2cpm.ord <- log2cpm.all[,ord]
counts.ord <- counts[,ord]

theta.ord <- theta[ord]

indices.train <- partitions[[1]]

out.methods <- mclapply(1:nrow(log2cpm.ord), function(g) {
#out.methods <- mclapply(1:10, function(g) {

  yy <- log2cpm.ord[g,]
  yy.train <- yy[indices.train]
  xx.train <- theta.ord[indices.train]
  not.zero.train <- yy.train > 0

  yy.train.nonzero <- yy.train[not.zero.train]
  xx.train.nonzero <- xx.train[not.zero.train]

  # local linear kernel estimate
  fit.npcirc.ll <- kern.reg.circ.lin(x = xx.train.nonzero,
                                     y = yy.train.nonzero, method="LL")
  npll.lfun <- approxfun(y=fit.npcirc.ll$y,x=fit.npcirc.ll$x, rule=2)
  npll.pred.y.train <- npll.lfun(xx.train.nonzero)
  npll.mad.pred <- mad(yy.train.nonzero-npll.pred.y.train, constant = 1)
  npll.mad.constant <- mad(yy.train.nonzero-mean(yy.train.nonzero), constant = 1)
  npll.mad.ratio <- npll.mad.pred/npll.mad.constant
  npll.pve <- 1-var(yy.train.nonzero-npll.pred.y.train)/var(yy.train.nonzero)

  # NW kernel estimate
  fit.npcirc.nw <- kern.reg.circ.lin(x = xx.train.nonzero,
                                     y = yy.train.nonzero, method="NW")
  npnw.lfun <- approxfun(y=fit.npcirc.nw$y,x=fit.npcirc.nw$x, rule=2)
  npnw.pred.y.train <- npnw.lfun(xx.train.nonzero)
  npnw.mad.pred <- mad(yy.train.nonzero-npnw.pred.y.train, constant = 1)
  npnw.mad.constant <- mad(yy.train.nonzero-mean(yy.train.nonzero), constant = 1)
  npnw.mad.ratio <- npnw.mad.pred/npnw.mad.constant
  npnw.pve <- 1-var(yy.train.nonzero-npnw.pred.y.train)/var(yy.train.nonzero)

  # smash gaussian
  is.zeros <- which(yy.train == 0)
  yy.train.impute <- yy.train
  for (i in 1:length(is.zeros)) {
    pos <- is.zeros[i]
    impute.val <- mean(yy.train[-is.zeros])
    yy.train.impute[pos] <- impute.val
  }
  fit.smash <- smash(x = yy.train.impute, model="gaus")

  smash.mad.pred <- mad(yy.train.impute-fit.smash, constant = 1)
  smash.mad.constant <- mad(yy.train.impute-mean(yy.train.impute), constant = 1)
  smash.mad.ratio <- smash.mad.pred/smash.mad.constant

  # smash poisson
  yy.counts <- counts.ord[g,]
  yy.counts.train <- yy.counts[indices.train]
  xx.train <- theta.ord[indices.train]
  not.zero.train <- yy.counts.train > 0

  is.zeros <- which(yy.counts.train == 0)
  yy.counts.train.impute <- yy.counts.train
  for (i in 1:length(is.zeros)) {
    pos <- is.zeros[i]
    impute.val <- round(mean(yy.counts.train[-is.zeros]))
    yy.counts.train.impute[pos] <- impute.val
  }
  fit.smash.pois <- smash(x = yy.counts.train.impute, model="poiss")

  smash.pois.mad.pred <- mad(yy.counts.train.impute-fit.smash.pois, constant = 1)
  smash.pois.mad.constant <- mad(yy.counts.train.impute-mean(yy.counts.train.impute), constant = 1)
  smash.pois.mad.ratio <- smash.pois.mad.pred/smash.pois.mad.constant

  return(list(smash.mad.ratio=smash.mad.ratio,
              smash.pois.mad.ratio=smash.pois.mad.ratio,
              npll.mad.ratio=npll.mad.ratio,
              npnw.mad.ratio=npnw.mad.ratio,
              npll.pve=npll.pve,
              npnw.pve=npnw.pve,
              yy.train=yy[indices.train],
              xx.train=theta.ord[indices.train],
              smash.xx=theta.ord[indices.train],
              smash.yy=fit.smash,
              smash.yy.pois=fit.smash.pois,
              npll.xx=xx.train.nonzero,
              npll.yy=npll.pred.y.train,
              npnw.xx=xx.train.nonzero,
              npnw.yy=npnw.pred.y.train))
}, mc.cores=10)

saveRDS(out.methods, file = file.path(dir,"output/npreg.Rmd/out.methods.rds"))


