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


library(Biobase)
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
theta <- pdata$theta

# select genes with few missings
ii.detected <- rowMeans(counts > 0) > .8
zz <- log2cpm.all[ii.detected, order(pdata$theta)]
pdata.ord <- pdata[order(pdata$theta),]

# source codes
source(file.path(dir,"code/utility.R"))
source(file.path(dir,"code/npreg/npreg.methods.R"))

chip_ids <- c("NA18511","NA18855","NA18870","NA19098","NA19101","NA19160")
out.methods <- vector("list", length(chip_ids))
names(out.methods) <- chip_ids
G <- nrow(zz)


## fit for training set and testing set
for (i in 1:length(out.methods)) {

  out.methods[[i]] <- mclapply(1:G, function(g) {
#  out.methods[[i]] <- mclapply(1:10, function(g) {

    yy <- zz[g,]
    out.trend <- fit.trendfilter(yy, c(1:length(yy)))

    for (perm in 1:nperm) {
      yy.perm <- sample(yy)
      out.trend <- fit.trendfilter(yy, c(1:length(yy)))
      f <- approxfun(x=out.trend$trend.pos, y=out.trend$trend.yy, rule=2)
      test.mad.pred <- mad(yy.test[nonzero.test]-f(pos.yy.test[nonzero.test]), constant = 1)
      test.mad.constant <- mad(yy.test[nonzero.test]-mean(yy.test[nonzero.test]), constant = 1)
      test.mad.ratio <- test.mad.pred/test.mad.constant

    }




    return(list(yy=yy,
                yy.train=yy.train,
                pos.yy.train=pos.yy.train,
                yy.test=yy.test,
                pos.test=pos.yy.test,
                pred.train.yy=out.trend$trend.yy,
                pred.train.xx=out.trend$trend.pos,
                pred.test.yy=f(pos.yy.test),
                pred.test.xx=pos.yy.test,
                test.mad.ratio=test.mad.ratio,
                train.mad.ratio=out.trend$trend.mad.ratio))
  }, mc.cores=5)
}

saveRDS(out.methods, file = file.path(dir,"output/npreg-trendfilter.Rmd/out.methods.rds"))


### fit for training set and testing set
# for (i in 1:length(out.methods)) {
#
#   id <- which(pdata.ord$chip_id == chip_ids[i])
#
#   out.methods[[i]] <- mclapply(1:G, function(g) {
# #  out.methods[[i]] <- mclapply(1:10, function(g) {
#
#     yy <- zz[g,]
#     yy.train <- yy[-id]
#     yy.test <- yy[id]
#     pos.yy.train <- c(1:length(yy))[-id]
#     pos.yy.test <- c(1:length(yy))[id]
#     nonzero.test <- yy.test!=0
#
#     out.trend <- fit.trendfilter(yy.train, pos.yy.train)
#
#     f <- approxfun(x=out.trend$trend.pos, y=out.trend$trend.yy, rule=2)
#
#     test.mad.pred <- mad(yy.test[nonzero.test]-f(pos.yy.test[nonzero.test]), constant = 1)
#     test.mad.constant <- mad(yy.test[nonzero.test]-mean(yy.test[nonzero.test]), constant = 1)
#     test.mad.ratio <- test.mad.pred/test.mad.constant
#
#     return(list(yy=yy,
#                 yy.train=yy.train,
#                 pos.yy.train=pos.yy.train,
#                 yy.test=yy.test,
#                 pos.test=pos.yy.test,
#                 pred.train.yy=out.trend$trend.yy,
#                 pred.train.xx=out.trend$trend.pos,
#                 pred.test.yy=f(pos.yy.test),
#                 pred.test.xx=pos.yy.test,
#                 test.mad.ratio=test.mad.ratio,
#                 train.mad.ratio=out.trend$trend.mad.ratio))
#   }, mc.cores=5)
#
# }
#
# saveRDS(out.methods, file = file.path(dir,"output/npreg-trendfilter.Rmd/out.methods.rds"))


