library(smashr)

set.seed(91)
ii.sub <- sample(length(Y_cdk1), 512, replace=F)
fit_smash <- smash(Y_cdk1[ii.sub])

library(splines)
fit_bspline <- lm(Y_cdk1 ~ bs(x = theta_train, df = 5))
pred_bspline <- predict(fit_bspline)

par(mfrow=c(1,1))
plot(pred_bspline)
plot(Y_cdk1)


library(CountClust)
df <- readRDS(file="../data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), which(pdata$chip_id == "NA18855")]
pdata <- pdata[which(pdata$chip_id == "NA18855"),]

fitgom <- FitGoM(t(counts[which(rownames(counts) %in% rownames(sig.genes)),]), K = 2:8)
sapply(fitgom, "[[", "BF")
all.equal(rownames(pdata), rownames(fitgom$clust_8$omega))
par(mfrow=c(2,4))
for (i in 1:8) {
  plot(pdata$theta, fitgom$clust_8$omega[,i])
}

plot(fitgom$clust_8$omega[,1],
     fitgom$clust_8$omega[,2])

df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)
counts <- exprs(df)[grep("ENSG", rownames(df)),]

Y=t(log2(1+(10^6)*(t(counts)/pdata$molecules)))
Y=log2cpm.quant
Y.sig <- Y[which(rownames(Y) %in% rownames(sig.genes)),]

Y.sig.correct <- do.call(rbind, lapply(1:nrow(Y.sig), function(g) {
  residuals <- resid(lm(Y.sig[g,] ~ factor(pdata$chip_id)))
  return(residuals)
}))
rownames(Y.sig.correct) <- rownames(Y.sig)

par(mfrow=c(2,3))
for (id in unique(pdata$chip_id)) {
  Y.sub <- Y.sig[,which(pdata$chip_id==id)]
  pr <- prcomp(t(Y.sub))
  plot(pr$x[,1], pr$x[,2])
}

par(mfrow=c(1,1))
plot(pdata$gfp.median.log10sum.adjust,
     pdata$rfp.median.log10sum.adjust)



install.packages("homals")
library(homals)
?homals

# writing new version of insample prediction
#' @title Estimate cell cycle ordering in the current sample
#'
#' @export
cycle.npreg.insample <- function(Y, theta, nbins=NULL, ncores=12,...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial <- cycle.npreg.mstep(Y = Y_ordered,
                               theta = theta_ordered_initial, ncores = ncores)

  # compute log-likelihood of the final fit
  mstep <- cycle.npreg.loglik(Y = initial$Y,
                              theta = initial$theta,
                              mu_est=initial$mu_est,
                              sigma_est=initial$sigma_est)
  out <- list(Y_ordered=mstep$Y,
              cell_times_est=mstep$theta,
              loglik_est=mstep.loglik,
              mu_est=mstep$mu_est,
              sigma_est=mstep$sigma_est)
}
