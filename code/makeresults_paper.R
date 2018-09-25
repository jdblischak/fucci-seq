# RNA-seq data quality, PCA, etc ---------------------------------------------------------

eset <- readRDS("data/eset-final.rds")

counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

# selection of technical factor
library(dplyr)
covariates <- pData(eset) %>% dplyr::select(experiment, well, chip_id,
                                                   concentration, raw:unmapped,
                                                   starts_with("detect"),  molecules)
# look at the first 6 PCs
pca_log2cpm <- prcomp(t(log2cpm.all), scale. = TRUE, center = TRUE)
pcs <- pca_log2cpm$x[, 1:6]


# R-square between PCs and the covariates
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

# plot heatmap
library(heatmap3)
heatmap3(r2, cexRow=1, cexCol=1, margins=c(2,12), scale="none",
         ylab="", main = "",
         Colv=F, showColDendro = F,
         labRow = c("C1: C1 batch",
                    "Well: C1 capture site (well) ID",
                    "Individual",
                    "cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"
                    ))

# r-squared between PC1 and individual
f1 <- lm(pcs[,1] ~ covariates$chip_id)
f2 <- lm(pcs[,1] ~ covariates$molecules)
summary(f1)

summary(f1)$adj.r.squared
summary(f2)$adj.r.squared


# PC proportions
100*((pca_log2cpm$sdev^2)/sum(pca_log2cpm$sdev^2))[1:6]


# Correlation between Technical factors
cor_tech <- cor(as.matrix(covariates[,4:11]),use="pairwise.complete.obs")
library(RColorBrewer)
heatmap3(cor_tech, symm = TRUE, margins=c(8,12),
         col=brewer.pal (9, "Blues" ), cexRow=1, cexCol=1, scale="none",
         labRow = c("cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"),
         labCol = c("cDNA",
                    "Reads raw",
                    "Reads UMI",
                    "Reads mapped",
                    "Reads unmapped",
                    "ERCC",
                    "ENSG",
                    "Molecules"))

summary(lm(covariates$detect_hs ~ covariates$experiment))

summary(lm(covariates$detect_hs ~ covariates$mapped))

summary(lm(covariates$detect_hs ~ covariates$chip_id))

summary(lm(covariates$detect_hs ~ covariates$molecules))

summary(lm(covariates$detect_hs ~ covariates$concentration))

# pca of top 10% experssing genes
log2cpm_mean <- rowMeans(log2cpm.all)
log2cpm_top <- log2cpm.all[rank(log2cpm_mean) / length(log2cpm_mean) > 1 - 0.1, ]
dim(log2cpm_top)

pca_top <- prcomp(t(log2cpm_top), scale. = T, center = T)

## look at the first 6 PCs
pcs <- pca_top$x[, 1:6]

## generate the data
r2_top <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
                 dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2_top[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

## plot heatmap
heatmap3(r2_top, cexRow=1, cexCol=1, margins=c(8,8),
         Colv=F, showColDendro = F,
         labRow = c("C1: C1 batch",
                    "Well: C1 capture site (well) ID",
                    "Individual",
                    "cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"
         ))



# Intensity batch effect ---------------------------------------------------------
lm.rfp <- lm(rfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.gfp <- lm(gfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.dapi <- lm(dapi.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
              data = pdata)

library(ibd)
aov.lm.rfp <- Anova(lm.rfp, type = "III")
aov.lm.gfp <- Anova(lm.gfp, type = "III")
aov.lm.dapi <- Anova(lm.dapi, type = "III")

aov.lm.rfp
aov.lm.gfp

library(ggplot2)
ggplot(pdata, aes(gfp.median.log10sum, col=factor(experiment))) +
  geom_density()

  facet_wrap(~ factor(chip_id))

ggplot(pdata, aes(x=factor(experiment), y=gfp.median.log10sum.adjust)) +
  geom_boxplot() + facet_wrap(~ factor(chip_id))

ggplot(pdata, aes(x=factor(experiment), y=gfp.median.log10sum.adjust)) +
  geom_boxplot()


# Fucci phase making and properties --------------------------------------------------
par(mfrow=c(1,1))
plot(x=pData(df)$rfp.median.log10sum.adjust, xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),
     y=pData(df)$gfp.median.log10sum.adjust, pch=16, cex=.5, col="gray50",
     xlab="RFP log10 sum intensity",
     ylab="GFP log10 sum intensity", axes=F)
axis(1); axis(2)

df <- readRDS("data/eset-final.rds")
pca <- prcomp(cbind(pData(df)$rfp.median.log10sum.adjust,
                    pData(df)$gfp.median.log10sum.adjust))
(pca$sdev^2)/sum(pca$sdev^2)
plot(pca$x[,1], pca$x[,2], pch=16, cex=.5, xlim=c(-4, 4), ylim=c(-4,4),
     xlab="PC1 (67%)", ylab="PC2 (33%)",
     main = "fucci intensities PC1 vs PC2", col="gray50", axes=F)
axis(1);axis(2)
abline(h=0,v=0, col="gray50", lty=2)
par(new=TRUE)

theta <- coord2rad(pca$x)
library(circular)
plot(circular(theta), stack=T, shrink=1.3, cex=.5, bins=200)


# fucci intensity and expression by phase
fits_all <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
genes <- names(pve_all_ord)[1:5]
data <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")

source("peco/R/fit.trendfilter.generic.R")
data_quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
theta_final <- (2*pi-as.numeric(theta))%%(2*pi)
sample_ord <- rownames(pData(df))[order(theta_final)]
data_quant_ord <- data_quant[,match(sample_ord,colnames(data_quant))]
fits_tmp <- lapply(1:5, function(i) {
  ii <- which(rownames(data_quant_ord)==genes[i])
  foo <- fit.trendfilter.generic(data_quant_ord[ii,])
  return(foo)
})
names(fits_tmp) <- genes

par(mfcol=c(6,1), mar=c(2,3,1,1))
plot(x=(2*pi-as.numeric(theta))%%(2*pi),
     y=pData(df)$gfp.median.log10sum.adjust, col="forestgreen",
     ylim=c(-1.5, 1.5), pch=16, cex=.5,
     xlab="Fucci phase", ylab="Fucci intensities adjusted for batch effect",
     main="Fucci intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi), labels=F)
abline(h=0, col="gray50", lty=2)
points(x=(2*pi-as.numeric(theta))%%(2*pi),
       y=pData(df)$rfp.median.log10sum.adjust, col="firebrick",
       ylim=c(-1.5, 1.5), pch=16, cex=.5)
for (i in 1:5) {
  plot(data_quant_ord[rownames(data_quant_ord)==genes[i],], col="gray50",
       xlab="Fucci phase",
       ylab="Normalized log2CPM", axes=F, cex=.7,
       main = labs[i])
  points(fits_tmp[[i]]$trend.yy, col="blue", pch=16, cex=.5)
  axis(2) #axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi), labels=F)
  abline(h=0, col="gray50", lty=2)
}


hist((2*pi-as.numeric(theta)), nclass=25,
     main="", xlab="Fucci phase")


# Top cyclical genes before training ----------------------------------------------

# import results
fits_all <- readRDS("../output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)

pve_all_ord <- pve_all[order(pve_all, decreasing = T)]

head(pve_all_ord)

genes <- names(pve_all_ord)[1:5]
data <- readRDS("../output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")
par(mfrow=c(2,3))
for (i in 1:5) {
  ii <- which(names(fits_all)==genes[i])
  plot(data[rownames(data)==genes[i],], col="gray50",
       xlab="Fucci ordering",
       ylab="Normalized log2CPM", axes=F,
       main = labs[i])
  axis(1); axis(2)
  points(fits_all[[ii]]$trend.yy, col="red", pch=16)
}


# Top cyclical genes in other datasets ----------------------------------------









# prediction error -------------------------------------------------------------
library(ggplot2)

double_diff_time_ind <- readRDS("output/method-train-summary-output.Rmd/double_diff_time_ind.rds")

double_diff_time_ind <- do.call(rbind, double_diff_time_ind)
double_diff_time_ind <- subset(double_diff_time_ind, methods=="supervised")

df <- rbind(data.frame(subset(double_diff_time_mixed, methods == "supervised"),
                       ind="mixed"),
            data.frame(subset(double_diff_time_mixed_permute, methods == "supervised"),
                       ind="mixed_permute"),
            double_diff_time_ind)

df$ind <- as.character(df$ind)
df$ind <- factor(df$ind, levels=c("mixed", "mixed_permute",
                                  "NA19098", "NA18855",
                                  "NA19101", "NA18511", "NA18870", "NA19160"),
                 labels=c("Mixed", "Mixed permute", LETTERS[1:6]))

# df$type <- df$ind
# df$type[df$ind=="Mixed"]

ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=diff_mean, group=ind, alpha=.05)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set") +
  ylim(0,.25) +
  geom_line(data=subset(df, ngenes <= 700 & ind=="Mixed"),
            aes(x=ngenes, y=diff_mean)) +
  geom_hline(yintercept=.25, col="gray50") +
  theme_light()



ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=(2^diff_se), group=ind, alpha=.05)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Mean squared error") + xlab("Top X cyclical genes") +
  ggtitle("Standard error of predictions") +
  #  ylim(0,.25) +
  geom_line(data=subset(df, ngenes <= 700 & ind=="Mixed"),
            aes(x=ngenes, y=(2^diff_se))) +
  #  geom_hline(yintercept=.25, col="gray50") +
  theme_light()



ggplot(subset(df, ngenes <= 400),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set up to top 400 genes") +
  ylim(.10,.17) +
  theme_light()



ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_errorbar(subset(df, ngenes <= 700),
                mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set") +
  ylim(0,.25)

ggplot(subset(df, ngenes == 5),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("using top 5 cyclical genes")



library(ggplot2)

double_diff_time_ind <- do.call(rbind, double_diff_time_ind)
double_diff_time_ind <- subset(double_diff_time_ind, methods=="supervised")

df <- rbind(data.frame(subset(double_diff_time_mixed, methods == "supervised"),
                       ind="mixed"),
            data.frame(subset(double_diff_time_mixed_permute, methods == "supervised"),
                       ind="mixed_permute"),
            double_diff_time_ind)

df$ind <- as.character(df$ind)
df$ind <- factor(df$ind, levels=c("mixed", "mixed_permute",
                                  "NA19098", "NA18855",
                                  "NA19101", "NA18511", "NA18870", "NA19160"),
                 labels=c("Mixed", "Mixed permute", LETTERS[1:6]))

# df$type <- df$ind
# df$type[df$ind=="Mixed"]

ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=diff_mean, group=ind, alpha=.05)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set") +
  ylim(0,.25) +
  geom_line(data=subset(df, ngenes <= 700 & ind=="Mixed"),
            aes(x=ngenes, y=diff_mean)) +
  geom_hline(yintercept=.25, col="gray50") +
  theme_light()



ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=(2^diff_se), group=ind, alpha=.05)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Mean squared error") + xlab("Top X cyclical genes") +
  ggtitle("Standard error of predictions") +
  #  ylim(0,.25) +
  geom_line(data=subset(df, ngenes <= 700 & ind=="Mixed"),
            aes(x=ngenes, y=(2^diff_se))) +
  #  geom_hline(yintercept=.25, col="gray50") +
  theme_light()



ggplot(subset(df, ngenes <= 400),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  # geom_errorbar(subset(df, ngenes <= 700),
  #               mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set up to top 400 genes") +
  ylim(.10,.17) +
  theme_light()



ggplot(subset(df, ngenes <= 700),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_errorbar(subset(df, ngenes <= 700),
                mapping=aes(ymin=diff_mean-diff_se, ymax=diff_mean+diff_se, col=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("Prediction error in the training set") +
  ylim(0,.25)

ggplot(subset(df, ngenes == 5),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Error margin (% arc length)") + xlab("Top X cyclical genes") +
  ggtitle("using top 5 cyclical genes")

# new version
ggplot(subset(df, ngenes <=50 & ind != "Mixed"),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  ylim(.10, .18) + xlim(0, 55) +
  geom_line(aes(col=ind), lty=2, lwd=.3) +
  geom_point(aes(color=ind)) + #geom_line(lty=3) +
  ylab("Prediction error (percent circle)") + xlab("Top X cyclical genes") +
  theme_light()



############################################################




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



# genes in oscope?
which(oscope$hgnc == "CDK1")
which(oscope$hgnc == "TOP2A")
which(oscope$hgnc == "UBE2C")
which(oscope$hgnc == "HIST1H4C")
which(oscope$hgnc == "HIST1H4E")






