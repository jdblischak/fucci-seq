conditioned on cell cycle state, we can model expression with projected multivarite linear model

take one gene,
split to different states,
extract ML parameters
reconstruct the distribution




##### Univariate spml to recover peak time

pdata.adj <- readRDS("../output/images-normalize-anova.Rmd/pdata.adj.rds")

pc.fucci <- prcomp(subset(pdata.adj, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)


# plot CDK1 vs cell time
macosko[macosko$hgnc == "CDK1",]
macosko[macosko$hgnc == "CDC6",]
cdk1 <- log2cpm.all[rownames(log2cpm.all)=="ENSG00000170312",]
cdc6 <- log2cpm.all[rownames(log2cpm.all)=="ENSG00000094804",]

# remove zeros
not.zeros <- cdk1 > 0 & cdc6 > 0
cdk1 <- cdk1[not.zeros]
cdc6 <- cdc6[not.zeros]
Theta.fucci <- 2*pi - Theta.fucci
Theta.fucci <- Theta.fucci[not.zeros]
Theta.cart <- Theta.cart[not.zeros,]

# fit spherical regression
library(Rfast)
est <- function(fit, cov) {
  mn.pt1 <- atan((fit$be[1,2] + fit$be[2,2]*cov)/(fit$be[1,1] + fit$be[2,1]*cov)) 
  mn.pt2 <- pi*(fit$be[1,1] + fit$be[2,1]*cov < 0)
  return((mn.pt1 + mn.pt2)%%(2*pi))
}


fit.cdk1 <- spml.reg(y=as.numeric(Theta.fucci), 
                  x=cdk1, seb = TRUE)
est.cdk1 <- est(fit.cdk1, cdk1)

fit.cdc6 <- spml.reg(y=as.numeric(Theta.fucci), 
                     x=cdc6, seb = TRUE)
est.cdc6 <- est(fit.cdc6, cdc6)

par(mfrow=c(2,2))
plot(as.numeric(Theta.fucci), cdk1)
plot(est.cdk1, cdk1)
plot(as.numeric(Theta.fucci), cdc6)
plot(est.cdc6, cdc6)


##### fit univariate spml to each gene

pdata.adj <- readRDS("../output/images-normalize-anova.Rmd/pdata.adj.rds")

pc.fucci <- prcomp(subset(pdata.adj, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)
Theta.fucci <- 2*pi - Theta.fucci

library(Rfast)
est <- function(fit, cov) {
  mn.pt1 <- atan((fit$be[1,2] + fit$be[2,2]*cov)/(fit$be[1,1] + fit$be[2,1]*cov)) 
  mn.pt2 <- pi*(fit$be[1,1] + fit$be[2,1]*cov < 0)
  return((mn.pt1 + mn.pt2)%%(2*pi))
}

fit.spml <- vector("list", nrow(log2cpm.all))
names(fit.spml) <- rownames(log2cpm.all)
for (g in 1:nrow(log2cpm.all)) {
  yy <- log2cpm.all[g,]
  not.zeros <- yy > 0
  yy.nonzero <- yy[not.zeros]
  Theta.fucci.nonzero <- Theta.fucci[not.zeros]

  fit <- spml.reg(y=as.numeric(Theta.fucci.nonzero), 
                     x=yy.nonzero, seb = TRUE)
  est.out <- est(fit, yy.nonzero)

  fit.spml[[g]] <- data.frame(yy=yy.nonzero,
                              est=est.out)
}
peaks <- sapply(fit.spml, function(x) x[which.max(x[,1]),2])
names(peaks) <- rownames(log2cpm.all)  


par(mfrow=c(1,1))
plot(circular(peaks), stack = TRUE)

peaks.cart <- cbind(sin(peaks), cos(peaks))


# cluster cell time
library(movMF)
peaks.res <- lapply(2:5, function(k) {
  movMF(peaks.cart, k=k, nruns = 100, kappa = list(common = TRUE))
})
k.list <- sapply(peaks.res, function(x) length(x$theta) + length(x$alpha) + 1)
bic <- sapply(1:length(peaks.res), function(i) {
  x <- peaks.res[[i]]
  k <- k.list[i]
  n <- nrow(peaks.cart)
  -2*x$L + k *(log(n) - log(2*pi)) })
plot(bic)
labs <- predict(peaks.res[[2]])


plot(circular(peaks), stack=TRUE)
points(circular(peaks)[labs==1], col = "red")
points(circular(peaks)[labs==2], col = "forestgreen")
points(circular(peaks)[labs==3], col = "blue")





##### fit univariate spml to genes jointly

pdata.adj <- readRDS("../output/images-normalize-anova.Rmd/pdata.adj.rds")

pc.fucci <- prcomp(subset(pdata.adj, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)
Theta.fucci <- 2*pi - Theta.fucci

library(Rfast)
est <- function(fit, cov) {
  mn.pt1 <- atan((fit$be[1,2] + fit$be[2,2]*cov)/(fit$be[1,1] + fit$be[2,1]*cov)) 
  mn.pt2 <- pi*(fit$be[1,1] + fit$be[2,1]*cov < 0)
  return((mn.pt1 + mn.pt2)%%(2*pi))
}

log2cpm.detected <- log2cpm.all[rowMeans(log2cpm.all>0)>.98,]

# impute missing values
log2cpm.all.impute <- matrix(0, 
                             nrow=nrow(log2cpm.detected), 
                             ncol=ncol(log2cpm.detected))
for (g in 1:nrow(log2cpm.detected)) {
  log2cpm.all.impute[g,] <- log2cpm.detected[g,]
  is.zero <- log2cpm.all.impute[g,]==0
  log2cpm.all.impute[g,is.zero] <- mean(log2cpm.detected[g,-which(is.zero)])
}

fit.spml <- spml.reg(y=as.numeric(Theta.fucci), 
                  x=t(log2cpm.all.impute), se=TRUE)

library(ashr)
beta1 <- ash(fit.spml$be[,1],fit.spml$seb[,1])
beta2 <- ash(fit.spml$be[,2],fit.spml$seb[,2])

summary(beta1$result$svalue)
summary(beta2$result$svalue)

hist(beta1$result$svalue, nclass=50)
hist(beta2$result$svalue, nclass=50)
plot(beta1$result$svalue, beta2$result$svalue)

table(beta1$result$svalue < .01, beta2$result$svalue < .01)
ii.sig <- which(beta1$result$svalue < .01 | beta2$result$svalue < .01)

plot(log2cpm.detected[ii.sig[1], order(Theta.fucci)])
plot(log2cpm.detected[ii.sig[2], order(Theta.fucci)])

sum(rownames(log2cpm.detected)[ii.sig] %in% macosko$ensembl)

tt <- (rownames(log2cpm.detected)[ii.sig])[rownames(log2cpm.detected)[ii.sig] %in% macosko$ensembl]

macosko[macosko$ensembl %in% tt,]

macosko[macosko$hgnc == "CDK1",]

# write.table(rownames(log2cpm.detected)[ii.sig],
#             file = "output/npreg.Rmd/spml.txt",
#             quote = F, row.names = F, col.names = F)









##### Univariate spml

pdata.adj <- readRDS("../output/images-normalize-anova.Rmd/pdata.adj.rds")

pc.fucci <- prcomp(subset(pdata.adj, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)


# plot CDK1 vs cell time
macosko[macosko$hgnc == "CDK1",]
cdk1 <- log2cpm.all[rownames(log2cpm.all)=="ENSG00000170312",]

# remove zeros
not.zeros <- cdk1 > 0
cdk1 <- cdk1[not.zeros]
Theta.fucci <- 2*pi - Theta.fucci
Theta.fucci <- Theta.fucci[not.zeros]
Theta.cart <- Theta.cart[not.zeros,]

# cluster cell time
library(movMF)
clust.res <- lapply(2:5, function(k) {
  movMF(Theta.cart, k=k, nruns = 100, kappa = list(common = TRUE))
})
k.list <- sapply(clust.res, function(x) length(x$theta) + length(x$alpha) + 1)
bic <- sapply(1:length(clust.res), function(i) {
  x <- clust.res[[i]]
  k <- k.list[i]
  n <- nrow(Theta.cart)
  -2*x$L + k *(log(n) - log(2*pi)) })
plot(bic)
labs <- predict(clust.res[[2]])

plot(Theta.fucci, stack=TRUE)
points(Theta.fucci[labs==1], col = "red")
points(Theta.fucci[labs==2], col = "forestgreen")
points(Theta.fucci[labs==3], col = "blue")


plot(y=cdk1, x=as.numeric(Theta.fucci), cex=1)
points(y=cdk1[labs==1], x=as.numeric(Theta.fucci)[labs==1], pch=16, cex=.7, col = "red")
points(y=cdk1[labs==2], x=as.numeric(Theta.fucci)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(y=cdk1[labs==3], x=as.numeric(Theta.fucci)[labs==3], pch=16, cex=.7, col = "blue")

summary(as.numeric(Theta.fucci)[labs==1])
summary(as.numeric(Theta.fucci)[labs==2])
summary(as.numeric(Theta.fucci)[labs==3])

# move the origin to 1.27
Theta.fucci.new <- vector("numeric", length(Theta.fucci))
cutoff <- min(Theta.fucci[labs==2])
Theta.fucci.new[Theta.fucci>=cutoff] <- Theta.fucci[Theta.fucci>=cutoff] - cutoff
Theta.fucci.new[Theta.fucci<cutoff] <- Theta.fucci[Theta.fucci<cutoff] - cutoff + 2*pi

plot(y=cdk1, x=as.numeric(Theta.fucci.new), cex=1)
points(y=cdk1[labs==1], x=as.numeric(Theta.fucci.new)[labs==1], pch=16, cex=.7, col = "red")
points(y=cdk1[labs==2], x=as.numeric(Theta.fucci.new)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(y=cdk1[labs==3], x=as.numeric(Theta.fucci.new)[labs==3], pch=16, cex=.7, col = "blue")

# fit spherical regression
library(Rfast)
est <- function(fit, cov) {
  mn.pt1 <- atan((fit$be[1,2] + fit$be[2,2]*cov)/(fit$be[1,1] + fit$be[2,1]*cov)) 
  mn.pt2 <- pi*(fit$be[1,1] + fit$be[2,1]*cov < 0)
  return((mn.pt1 + mn.pt2)%%(2*pi))
}

fit.1 <- spml.reg(y=as.numeric(Theta.fucci.new)[labs==1], 
         x=cdk1[labs==1], seb = TRUE)
fit.2 <- spml.reg(y=as.numeric(Theta.fucci.new)[labs==2], 
                  x=cdk1[labs==2], seb = TRUE)
fit.3 <- spml.reg(y=as.numeric(Theta.fucci.new)[labs==3], 
                  x=cdk1[labs==3], seb = TRUE)

est.1 <- est(fit.1, cdk1[labs==1])
est.2 <- est(fit.2, cdk1[labs==2])
est.3 <- est(fit.3, cdk1[labs==3])

par(mfrow=c(1,2))
plot(x=cdk1, y=as.numeric(Theta.fucci.new), cex=1)
points(x=cdk1[labs==1], y=as.numeric(Theta.fucci.new)[labs==1], pch=16, cex=.7, col = "red")
points(x=cdk1[labs==2], y=as.numeric(Theta.fucci.new)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(x=cdk1[labs==3], y=as.numeric(Theta.fucci.new)[labs==3], pch=16, cex=.7, col = "blue")
plot(x=cdk1[labs==1], y=est.1, cex=1, ylim = c(0,2*pi), xlim=range(cdk1), col = "red")
points(x=cdk1[labs==2], y=est.2, cex=1, col = "forestgreen")
points(x=cdk1[labs==3], y=est.3, cex=1, col = "blue")

### not conditioned on state
fit <- spml.reg(y=as.numeric(Theta.fucci.new), 
                  x=cdk1, seb = TRUE)
est.out <- est(fit, cdk1)
par(mfrow=c(1,2))
plot(x=cdk1, y=as.numeric(Theta.fucci.new), cex=1)
points(x=cdk1[labs==1], y=as.numeric(Theta.fucci.new)[labs==1], pch=16, cex=.7, col = "red")
points(x=cdk1[labs==2], y=as.numeric(Theta.fucci.new)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(x=cdk1[labs==3], y=as.numeric(Theta.fucci.new)[labs==3], pch=16, cex=.7, col = "blue")
plot(x=cdk1, y=est.out, cex=1, ylim = c(0,2*pi), xlim=range(cdk1))


### two covariates
pdata.adj <- readRDS("../output/images-normalize-anova.Rmd/pdata.adj.rds")

pc.fucci <- prcomp(subset(pdata.adj, 
                          select=c("rfp.median.log10sum.adjust",
                                   "gfp.median.log10sum.adjust")),
                   center = T, scale. = T)
Theta.cart <- pc.fucci$x
library(circular)
Theta.fucci <- coord2rad(Theta.cart)

macosko[macosko$hgnc == "CDK1",]
macosko[macosko$hgnc == "TOP2A",]
cdk1 <- log2cpm.all[rownames(log2cpm.all)=="ENSG00000170312",]
top2a <- log2cpm.all[rownames(log2cpm.all)=="ENSG00000131747",]

# remove zeros
not.zeros <- cdk1 > 0 & top2a > 0
cdk1 <- cdk1[not.zeros]
top2a <- top2a[not.zeros]
Theta.fucci <- 2*pi - Theta.fucci
Theta.fucci <- Theta.fucci[not.zeros]
Theta.cart <- Theta.cart[not.zeros,]

par(mfrow=c(1,2))
plot(y=cdk1, x=as.numeric(Theta.fucci.new), cex=1, main = "CDK1")
points(y=cdk1[labs==1], x=as.numeric(Theta.fucci.new)[labs==1], pch=16, cex=.7, col = "red")
points(y=cdk1[labs==2], x=as.numeric(Theta.fucci.new)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(y=cdk1[labs==3], x=as.numeric(Theta.fucci.new)[labs==3], pch=16, cex=.7, col = "blue")
plot(y=top2a, x=as.numeric(Theta.fucci.new), cex=1, main = "TOP2A")
points(y=top2a[labs==1], x=as.numeric(Theta.fucci.new)[labs==1], pch=16, cex=.7, col = "red")
points(y=top2a[labs==2], x=as.numeric(Theta.fucci.new)[labs==2], pch=16, cex=.7, col = "forestgreen")
points(y=top2a[labs==3], x=as.numeric(Theta.fucci.new)[labs==3], pch=16, cex=.7, col = "blue")


# fit spherical regression
library(Rfast)
est <- function(fit, cov) {
  mn.pt1 <- atan((fit$be[1,2] + fit$be[2,2]*cov)/(fit$be[1,1] + fit$be[2,1]*cov)) 
  mn.pt2 <- pi*(fit$be[1,1] + fit$be[2,1]*cov < 0)
  return((mn.pt1 + mn.pt2)%%(2*pi))
}

fit.cdk1 <- spml.reg(y=as.numeric(Theta.fucci.new), 
                  x=cdk1, seb = TRUE)
fit.top2a <- spml.reg(y=as.numeric(Theta.fucci.new), 
                  x=top2a, seb = TRUE)
est.cdk1 <- est(fit.cdk1, cdk1)
est.top2a <- est(fit.top2a, top2a)

par(mfrow=c(1,2))
plot(x=cdk1, y=est.cdk1, cex=1)
plot(x=top2a, y=est.top2a, cex=1)


fit.joint <- spml.reg(y=as.numeric(Theta.fucci.new), 
                   x=cbind(cdk1,top2a), seb = TRUE)


###############

par(mfrow=c(1,2))
plot(theta.sub, mol.sub)
gfp.sub <- pdata$gfp.median.log10sum[ind]
rfp.sub <- pdata$rfp.median.log10sum[ind]
dapi.sub <- pdata$dapi.median.log10sum[ind]
plot(-rev(theta.sub), rev(gfp.sub), col = "green")
points(-rev(theta.sub), rev(rfp.sub), col = "red")
points(-rev(theta.sub), rev(dapi.sub), col = "blue")