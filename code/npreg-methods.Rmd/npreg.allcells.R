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


## Load data
dir <- "/project2/gilad/joycehsiao/fucci-seq"
df <- readRDS(file.path(dir,"data/eset-filtered.rds"))
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- readRDS(file.path(dir,"output/seqdata-batch-correction.Rmd/log2cpm.rds"))
# log2cpm.adjust <- readRDS("../output/seqdata-batch-correction.Rmd/log2cpm.adjust.rds")
log2cpm <- log2cpm.all[grep("ENSG", rownames(log2cpm.all)), ]

# import corrected intensities
pdata.adj <- readRDS(file.path(dir,"output/images-normalize-anova.Rmd/pdata.adj.rds"))

macosko <- readRDS(file.path(dir,"data/cellcycle-genes-previous-studies/rds/macosko-2015.rds"))

log2cpm.detected <- log2cpm[rowMeans(log2cpm>0)>.8,]

#log2cpm.high <- log2cpm.detected[order(rowMeans(log2cpm.detected))[1:100], ]

source(file.path(dir,"code/utility.R"))


## Projected normal on PCs of Red/Green
pc.fucci <- prcomp(subset(pdata.adj,
                        select=c("rfp.median.log10sum.adjust",
                                 "gfp.median.log10sum.adjust")),
                 center = T, scale. = T)
Theta.cart <- pc.fucci$x
Theta.fucci <- coord2rad(Theta.cart)



## Compare different methods
#log2cpm.detected <- log2cpm[rowMeans(log2cpm>0)>.8,]

## order cells based on estimated time
# partitions <- partitionSamples.n(c(1:length(Theta.fucci)), runs=1,
#                                  c(512,400))[[1]][[1]]
ord <- order(as.numeric(Theta.fucci))
#log2cpm.detected.ord <- log2cpm.detected[,ord]
log2cpm.ord <- log2cpm[,ord]
theta.ord <- as.numeric(Theta.fucci)[ord]

#indices.train <- partitions[[1]]

out.methods <- mclapply(1:nrow(log2cpm.ord), function(g) {
#out.methods <- mclapply(1:10, function(g) {

  yy <- log2cpm.ord[g,]
  # yy.train <- yy[indices.train]
  # xx.train <- theta.ord[indices.train]
  yy.train <- yy
  xx.train <- theta.ord
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

  # # smash
  # is.zeros <- which(yy.train == 0)
  # yy.train.impute <- yy.train
  # for (i in 1:length(is.zeros)) {
  #   pos <- is.zeros[i]
  #   impute.val <- mean(yy.train[-is.zeros])
  #   yy.train.impute[pos] <- impute.val
  # }
  # fit.smash <- smash(x = yy.train.impute, model="gaus")
  #
  # smash.mad.pred <- mad(yy.train.impute-fit.smash, constant = 1)
  # smash.mad.constant <- mad(yy.train.impute-mean(yy.train.impute), constant = 1)
  # smash.mad.ratio <- smash.mad.pred/smash.mad.constant

  return(list(npll.mad.ratio=npll.mad.ratio,
              npnw.mad.ratio=npnw.mad.ratio,
              npll.pve=npll.pve,
              npnw.pve=npnw.pve,
              yy.train=yy.train,
              xx.train=theta.ord,
              npll.xx=xx.train.nonzero,
              npll.yy=npll.pred.y.train,
              npnw.xx=xx.train.nonzero,
              npnw.yy=npnw.pred.y.train))
}, mc.cores=10)

saveRDS(out.methods, file = file.path(dir,"output/npreg.Rmd/out.methods.allcells.rds"))


