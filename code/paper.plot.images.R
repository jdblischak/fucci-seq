df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

macosko <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")

theta <- 2*pi - pdata$theta
log2cpm.all <- log2cpm.all[,order(theta)]
pdata <- pdata[order(theta),]
theta <- theta[order(theta)]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")


red_high <- rownames(pdata)[which.max(pdata$rfp.median.log10sum.adjust)]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        rfp.median.log10sum.adjust)), theta)[which(rownames(pdata) %in% red_high),]

red_low <- rownames(pdata)[which.min(pdata$rfp.median.log10sum.adjust)]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        rfp.median.log10sum.adjust)), theta)[which(rownames(pdata) %in% red_low),]

red_middle <- rownames(pdata)[pdata$rfp.median.log10sum.adjust < .6 & pdata$rfp.median.log10sum.adjust > .5 & theta > 3 & theta < 3.01]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        rfp.median.log10sum.adjust)), theta)[which(rownames(pdata) %in% red_middle),]

green_high <- rownames(pdata)[which.max(pdata$gfp.median.log10sum.adjust)]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        gfp.median.log10sum.adjust)), theta)[which(rownames(pdata) %in% green_high),]

# green_low <- rownames(pdata)[which.min(pdata$gfp.median.log10sum.adjust)]
# cbind(with(pdata, cbind(rownames(pdata),
#                         image_individual, image_label,
#                         gfp.median.log10sum.adjust)), theta)[which(rownames(pdata) %in% green_low),]

theta_min <- rownames(pdata)[which.min(theta)]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        gfp.median.log10sum.adjust)), theta)[which.min(theta),]

theta_max <- rownames(pdata)[which.max(theta)]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        gfp.median.log10sum.adjust)), theta)[which.max(theta),]

others1 <- rownames(pdata)[which(rownames(pdata)=="20170906-E07")]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        rfp.median.log10sum.adjust,
                        gfp.median.log10sum.adjust)), theta)[which(rownames(pdata)=="20170906-E07"),]

#rownames(pdata)[which(theta < 3.5 & theta > 3.46)]
others2 <- rownames(pdata)[which(rownames(pdata)=="20170921-D05")]
cbind(with(pdata, cbind(rownames(pdata),
                        image_individual, image_label,
                        rfp.median.log10sum.adjust,
                        gfp.median.log10sum.adjust)), theta)[which(rownames(pdata)=="20170921-D05"),]



plot(x=theta, pdata$gfp.median.log10sum.adjust, col = "forestgreen", cex=.6, pch=16,
     ylim = c(-1.8,1.5),
     ylab = "log10 sum of adjusted pixel intensities",
     xlab = "Estimated cell time based on FUCCI intensities",
     axes=F)
axis(1); axis(2)
#points(x=theta, pdata$rfp.median.log10sum.adjust, col = "red", cex=.6, pch=16)
abline(h=0, col="gray50", lty=3)
points(x=theta[which(rownames(pdata) == red_low)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == red_low)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == red_high)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == red_high)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == red_middle)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == red_middle)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == green_high)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == green_high)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == theta_min)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == theta_min)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == theta_max)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == theta_max)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == others1)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == others1)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == others2)],
       y= pdata$gfp.median.log10sum.adjust[which(rownames(pdata) == others2)],
       pch=5, col="black", cex=1, lwd=1)


plot(x=theta, pdata$rfp.median.log10sum.adjust, col = "red", cex=.6, pch=16,
     ylim = c(-1.8,1.5),
     ylab = "log10 sum of adjusted pixel intensities",
     xlab = "Estimated cell time based on FUCCI intensities",
     axes=F)
axis(1); axis(2)
#points(x=theta, pdata$rfp.median.log10sum.adjust, col = "red", cex=.6, pch=16)
abline(h=0, col="gray50", lty=3)
points(x=theta[which(rownames(pdata) == red_low)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == red_low)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == red_high)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == red_high)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == red_middle)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == red_middle)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == green_high)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == green_high)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == theta_min)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == theta_min)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == theta_max)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == theta_max)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == others1)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == others1)],
       pch=5, col="black", cex=1, lwd=1)
points(x=theta[which(rownames(pdata) == others2)],
       y= pdata$rfp.median.log10sum.adjust[which(rownames(pdata) == others2)],
       pch=5, col="black", cex=1, lwd=1)


cbind(rownames(pdata),
      pdata$image_individual,
      pdata$image_label, theta)[which(rownames(pdata) %in% c(red_high,
                                                             red_middle, red_low,
                                                             green_high,
                                                             theta_min, theta_max,
                                                             others1,
                                                             others2)),]

# cdk1 <- with(macosko, ensembl[which(hgnc=="CDK1")])
# plot(log2cpm.quant[rownames(log2cpm.quant)==cdk1,])
#
# cdc6 <- with(macosko, ensembl[which(hgnc=="CDC6")])
# plot(log2cpm.quant[rownames(log2cpm.quant)==cdc6,])
