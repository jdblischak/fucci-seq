Y <- Y_train
library(circular)
pc_res <- prcomp(t(Y), scale = TRUE)

plot(pc_res$x[,1],
     pc_res$x[,2])

theta_train_pca <- initialize_cell_times(Y_train)

cdk1 <- macosko$ensembl[which(macosko$hgnc=="CDK1")]

Y_cdk1 <- Y[rownames(Y) == cdk1, ]

par(mfrow=c(2,1))
plot(theta_train_pca, Y_cdk1[match(names(theta_train_pca),names(Y_cdk1))])
plot(theta_train, Y_cdk1[match(names(theta_train),names(Y_cdk1))])

par(mfrow=c(1,1))
plot(theta_train, theta_train_pca)

all.equal(names(theta_train_pca), rownames(pdata[part_indices[[1]]$train,]))
ii.match <- match(rownames(pdata[part_indices[[1]]$train,]),
                  names(theta_train_pca))

par(mfcol=c(3,2))
plot(theta_train_pca[ii.match], pdata[part_indices[[1]]$train,]$rfp.median.log10sum.adjust,
     col="red")
plot(theta_train_pca[ii.match], pdata[part_indices[[1]]$train,]$gfp.median.log10sum.adjust,
     col="green")
plot(theta_train_pca[ii.match], pdata[part_indices[[1]]$train,]$dapi.median.log10sum.adjust,
     col="blue", ylim = c(-1,1))

plot(theta_train[ii.match], pdata[part_indices[[1]]$train,]$rfp.median.log10sum.adjust,
     col="red")
plot(theta_train[ii.match], pdata[part_indices[[1]]$train,]$gfp.median.log10sum.adjust,
     col="green")
plot(theta_train[ii.match], pdata[part_indices[[1]]$train,]$dapi.median.log10sum.adjust,
     col="blue", ylim = c(-1,1))


head(names(theta_train_pca))
head(names(theta_train))

head((theta_train_pca))
head((theta_train))

plot(Y_cdk1[match(names(theta_train),names(Y_cdk1))])

plot(theta_train_pca, theta_train)
