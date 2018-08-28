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


# more recent figures ------------------------------------------------------

df <- readRDS("data/eset-final.rds")
pca <- prcomp(cbind(pData(df)$rfp.median.log10sum.adjust,
                    pData(df)$gfp.median.log10sum.adjust))
(pca$sdev^2)/sum(pca$sdev^2)
plot(pca$x[,1], pca$x[,2], pch=16, cex=.7, xlim=c(-1.2, 1.2), ylim=c(-1.2,1.2),
     xlab="PC1 (67%)", ylab="PC2 (33%)",
     main = "fucci intensities PC1 vs PC2")
abline(h=0,v=0, col="gray50", lty=2)
points(cos(pData(df)$theta), sin(pData(df)$theta), col="brown", pch=16, cex=.7)

plot(x=(pData(df)$theta -2)%%(2*pi),
     y=pData(df)$gfp.median.log10sum.adjust, col="forestgreen",
     ylim=c(-1.5, 1.5), pch=16, cex=.7,
     xlab="fucci time", ylab="fucci intensities adjusted for batch effect",
     main="Fucci intensities ordered by fucci time")
abline(h=0, col="gray50", lty=2)
points(x=(pData(df)$theta -2)%%(2*pi),
       y=pData(df)$rfp.median.log10sum.adjust, col="firebrick",
       ylim=c(-1.5, 1.5), pch=16, cex=.7)
hist((pData(df)$theta -2)%%(2*pi), nclass=25,
     main="fucci time frequency", xlab="fucci time")

theta <- (pData(df)$theta -2)%%(2*pi)
plot(pData(df)$theta,
     pData(df)$dapi.median.log10sum.adjust)


cdk1 <- rownames(fData(df))[which(fData(df)$name=="CDK1")]
dtl <- rownames(fData(df))[which(fData(df)$name=="DTL")]
expr_normed <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")

all.equal(colnames(expr_normed), rownames(pData(df)))
expr_normed2 <- expr_normed[,match(rownames(pData(df)),colnames(expr_normed))]

theta <- (pData(df)$theta -2)%%(2*pi)
cdc6 <- rownames(fData(df))[which(fData(df)$name=="CDC6")]
yy <- expr_normed2[rownames(expr_normed2)==cdc6,]
fit.cdc6 <- fit.trendfilter.generic(yy[order(theta)], polyorder=2)
plot(x=theta[order(theta)], y=yy[order(theta)],
     ylab="log2CPM normalized expression", xlab="fucci time",
     main = "CDC6", cex=.7, col="gray40")
abline(h=0, col="gray50")
points(x=theta[order(theta)],
       y=fit.cdc6$trend.yy, col ="brown", pch=16, cex=.6)


# theta <- (pData(df)$theta -2)%%(2*pi)
# yy <- expr_normed2[rownames(expr_normed2)==dtl,]
# fit.dtl <- fit.trendfilter.generic(yy[order(theta)], polyorder=2)
# plot(x=theta[order(theta)], y=yy[order(theta)],
#      ylab="log2CPM normalized expression", xlab="fucci time",
#      main = "DTL", cex=.7)
# points(x=theta[order(theta)],
#        y=fit.dtl$trend.yy, col ="brown", pch=16, cex=.6)

theta <- (pData(df)$theta -2)%%(2*pi)
yy <- expr_normed2[rownames(expr_normed2)==cdk1,]
fit.cdk1 <- fit.trendfilter.generic(yy[order(theta)], polyorder=2)
plot(x=theta[order(theta)], y=yy[order(theta)],
     ylab="log2CPM normalized expression", xlab="fucci time",
     main = "CDK1", cex=.7, col="gray40")
abline(h=0, col="gray50")
points(x=theta[order(theta)],
       y=fit.cdk1$trend.yy, col ="brown", pch=16, cex=.6)




# getting gene information
genes_list <- readRDS(file = "data/results/results_topgenes.rds")
genes_list_symbols <- readRDS("output/method-train-classifiers-genes.Rmd/genes_list_symbols.rds")


genes_list_symbols[[1]]
sapply(genes_list, length)[1:10]
sapply(genes_list_symbols, nrow)[1:10]


seurat.genes <- readLines(
  con = "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt")
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])

symbs <- c(genes_list_symbols[[1]]$hgnc_symbol[1:4], "HIST1H4E")
which(symbs %in% unlist(seurat.genes))
symbs[1:3]


symbs <- c(genes_list_symbols[[1]]$hgnc_symbol[1:4], "HIST1H4E")
which(symbs %in% unlist(seurat.genes))
symbs[ii]


symbs <- c(genes_list_symbols[[2]]$hgnc_symbol[which(!is.na(genes_list_symbols[[2]]$hgnc_symbol))], "HIST1H4E")
ii <- which(symbs %in% unlist(seurat.genes))
ii2 <- which(!(symbs %in% unlist(seurat.genes)))
symbs[ii]
symbs[ii2]


symbs <- c(genes_list_symbols[[3]]$hgnc_symbol[which(!is.na(genes_list_symbols[[3]]$hgnc_symbol))], "HIST1H4E")
ii <- which(symbs %in% unlist(seurat.genes))
ii2 <- which(!(symbs %in% unlist(seurat.genes)))
symbs[ii]
symbs[ii2]



symbs <- c(genes_list_symbols[[11]]$hgnc_symbol[which(!is.na(genes_list_symbols[[11]]$hgnc_symbol))],
           "HIST1H4E", "HIST1H4L", "HIST1H4B")

oo <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")

bb <- which(symbs %in% oo$hgnc | symbs %in% unlist(seurat.genes))
length(bb)
symbs[bb]
bb2 <- which(!(symbs %in% oo$hgnc | symbs %in% unlist(seurat.genes)))
length(bb2)
symbs[bb2]


symbs <- c(genes_list_symbols[[1]]$hgnc_symbol[1:4], "HIST1H4E")
bb <- which(symbs %in% oo$hgnc | symbs %in% unlist(seurat.genes))
length(bb)
bb2 <- which(!(symbs %in% oo$hgnc | symbs %in% unlist(seurat.genes)))
symbs[bb2]


