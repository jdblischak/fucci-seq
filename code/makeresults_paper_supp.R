# prediction error: training data specific top genes -----------------------------------

diff_time_wrapper <- function(results_list) {

  methods_list <- sapply(names(results_list),
                         function(x) strsplit(x, split=".", fixed=TRUE)[[1]][2])

  diff_time_list <- do.call(rbind, lapply(1:length(results_list), function(i) {
    diff_time <- results_list[[i]]$diff_time
    diff_mean <- mean(diff_time/2/pi)
    diff_mse <- sd(diff_time/2/pi)^2 #/sqrt(ncol(results_list[[1]]$Y))

    return(data.frame(diff_mean=diff_mean,
                      diff_mse=diff_mse,
                      methods=methods_list[i]))
  }) )
  return(diff_time_list)
}

# code saved in code/working/finalizing/job_run_methods.train.ind.R

ngenes <- c(5, seq(10,50, by=10))
#genes_list <- readRDS("data/results/finalizing/ind_NA18511_cyclical_genes.rds")
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")

eval_res <- do.call(rbind, lapply(1:length(inds), function(j) {
  ind <- inds[j]
  foo <- do.call(rbind, lapply(1:length(ngenes), function(i) {
    ngene <- ngenes[i]
    #  train_topX <- do.call(rbind, lapply(1:5, function(fold) {
    fl_name <- list.files("data/results/finalizing",
                          pattern=paste0("ind_",ind,"_results.top",ngene,".rds"),
                          full.names = TRUE)
    df <- readRDS(fl_name)
    out <- diff_time_wrapper(df$fit.test)
    out$ngenes <- ngene
    return(out)
  }) )
  foo$ind <- ind
  return(foo)
}) )


library(ggplot2)
ggplot(subset(eval_res,methods=="supervised"),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind), alpha=.5) + #geom_line(lty=3) +
  ylab("Prediction error (% circle)") + xlab("Top X cyclical genes") +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.25) +
  geom_hline(yintercept=.25, col="gray50") +
  theme_light()


ggplot(subset(eval_res,methods=="supervised"),
       aes(x=ngenes, y=diff_mse, group=ind)) +
  geom_point(aes(color=ind), alpha=.5) +
  ylim(0,.25) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylab("Mean squared error (% circle)") + xlab("Top X cyclical genes") +
  theme_light()


# mean prediction error betewen 5 to 100
eval_res_sub <- subset(eval_res,methods=="supervised")
aggregate(diff_mean ~ ngenes, data=eval_res_sub, mean)

# prediction error range of 5 genes predictor
subset(eval_res,methods=="supervised" & ngenes==5)


# check if the top 5 genes are the same across the folds (individual scenarios)
pvegenes <- lapply(1:length(inds), function(i) {
  ind <- inds[i]
  fl_name <- list.files("data/results/finalizing",
                        pattern=paste0("ind_",ind,"_cyclical_genes.rds"),
                        full.names = TRUE)
  df <- readRDS(fl_name)
  return(df)
})
names(pvegenes) <- inds

all.equal(rownames(pvegenes[[1]]), rownames(pvegenes[[2]]))
genes_list <- rownames(pvegenes[[1]])
topgenes <- lapply(1:length(pvegenes), function(i) {
  tmp <- genes_list[order(pvegenes[[i]]$pve, decreasing = T)]
  tmp[1:5]
})
names(topgenes) <- inds

# check if the top 5 in all individual scenarios are the same
fits_all <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
topgenes_all <- names(pve_all_ord)[1:5]

for (i in 1:length(pvegenes)) {
  print(sum(topgenes_all %in% topgenes[[i]]))
}

# check which gene is not in overall top 5
# KIF23: Kinase responsible for mitosis
# instead, HISTH4C is in the overall but not in NA19101 top 5
tops <- rownames(pvegenes[[4]])[order(pvegenes[[4]], decreasing=T)][1:5]
tops[!(tops %in% topgenes_all)]

pvegenes[[4]][rownames(pvegenes[[4]]) == "ENSG00000197061",]
head(pvegenes[[4]][order(pvegenes[[4]], decreasing=T),], 20)

pvegenes[[1]][rownames(pvegenes[[1]]) %in% topgenes[[i]],]
pvegenes[[2]][rownames(pvegenes[[2]]) %in% topgenes[[i]],]
pvegenes[[3]][rownames(pvegenes[[3]]) %in% topgenes[[i]],]
pvegenes[[4]][rownames(pvegenes[[4]]) %in% topgenes[[i]],]
pvegenes[[5]][rownames(pvegenes[[5]]) %in% topgenes[[i]],]
pvegenes[[6]][rownames(pvegenes[[6]]) %in% topgenes[[i]],]


# make figure 2
eval_res_sub2 <- subset(eval_res,
                        (ngenes <= 50) &( methods == "supervised"|methods=="seurat"))

eval_res_sub2 <- subset(eval_res_sub2,
                        !(methods == "seurat" & (ngenes > 5)))


eval_res_sub2$predictor <- NULL

eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==5] <- "peco 5 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==10] <- "peco 10 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==20] <- "peco 20 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==30] <- "peco 30 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==40] <- "peco 40 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==50] <- "peco 50 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="seurat" & eval_res_sub2$ngenes ==5] <- "seurat 97 genes"

eval_res_sub2$predictor <- factor(eval_res_sub2$predictor)
eval_res_sub2$predictor <- factor(eval_res_sub2$predictor,
                                  levels(eval_res_sub2$predictor)[c(5,1,2,3,4,6,7)])

ggplot(eval_res_sub2,
       aes(x=predictor, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind), alpha=.7) +
  stat_summary(fun.y=mean,geom="point",shape="_",stroke=10, group=1) +
  ylab("Predictor error (percent circle)") + xlab("Predictors") +
  theme_light() +
  theme(axis.text.x = element_text(angle=35, vjust=.5, hjust=.3))





# evaluate Seurat

labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")
genenames <- c("ENSG00000170312", "ENSG00000175063", "ENSG00000131747",
               "ENSG00000198518", "ENSG00000197061")
par(mfrow=c(1,5), mar=c(5,3,2,1))
# genenames <- rownames(eval_seurat[[1]]$`5`$fit.supervised$Y)
# genelabs <- c("CDK1", "UBE2C", "TOP2A", "H1")
for (g in 1:length(genenames)) {
  plot(x=eval_seurat[[1]]$`5`$fit.supervised$cell_times_est,
       y=eval_seurat[[1]]$`5`$fit.supervised$Y[rownames(fit.test$fit.supervised$Y)==genenames[g],],
       main=labs[g],
       ylab="Normalized log2CPM", pch="",
       xlab="Fucci phase")
  points(x=eval_seurat[[1]]$`5`$fit.supervised$cell_times_est,
         y=eval_seurat[[1]]$`5`$fit.supervised$Y[rownames(fit.test$fit.supervised$Y)==genenames[g],],
         col=c("darkgoldenrod1", "coral", "darkred")[eval_seurat[[1]]$`5`$fit.seurat$assignments],
         pch=16, cex=.7)
}



par(mfrow=c(1,5), mar=c(5,3,2,1))
# genenames <- rownames(eval_seurat[[1]]$`5`$fit.supervised$Y)
# genelabs <- c("CDK1", "UBE2C", "TOP2A", "H1")
for (g in 1:length(genenames)) {
  boxplot(eval_seurat[[1]]$`5`$fit.supervised$Y[rownames(fit.test$fit.supervised$Y)==genenames[g],]
          ~ eval_seurat[[1]]$`5`$fit.seurat$assignments,
          main=labs[g],
          col=c("darkgoldenrod1", "coral", "darkred"),
          ylab="Normalized log2CPM",
          xlab="Discrete class")
  # points(x=eval_seurat[[1]]$`5`$fit.supervised$cell_times_est,
  #        y=eval_seurat[[1]]$`5`$fit.supervised$Y[rownames(fit.test$fit.supervised$Y)==genenames[g],],
  #        col=c("darkgoldenrod1", "coral", "darkred")[eval_seurat[[1]]$`5`$fit.seurat$assignments],
  #        pch=16, cex=.7)
}

library(ggplot2)
df_tmp <- as.data.frame(do.call(cbind,eval_seurat[[1]]$`5`$fit.seurat))
df_tmp$assignments <- factor(df_tmp$assignments, levels=c("1", "2", "3"),
                             labels=c("G1", "S", "G2M"))
ggplot(df_tmp,
       aes(y=ref_time,x=factor(assignments),
           fill=assignments)) +
  geom_violin() +
  geom_boxplot(width=0.1) + theme_light() +
  scale_fill_manual(values=c("darkgoldenrod1", "coral", "darkred"))



# seurat continuous time
df_tmp <- as.data.frame(do.call(cbind,eval_seurat[[1]]$`5`$fit.seurat))
df_tmp$assignments <- factor(df_tmp$assignments, levels=c("1", "2", "3"),
                             labels=c("G1", "S", "G2M"))
ggplot(df_tmp,
       aes(y=pred_time_shift,x=factor(assignments),
           fill=assignments)) +
  ylab("Continuous phase based on expression scores") +
  xlab("Discrete cell-cycle class") +
  geom_violin() +
  geom_boxplot(width=0.1) + theme_light() +
  scale_fill_manual(values=c("darkgoldenrod1", "coral", "darkred"))


# Fucci phase by seurat classes results
plot(eval_seurat[[1]]$`5`$fit.seurat$ref_time,
     eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift,
     pch=16, cex=.7, xlab="")
library(circular)
par(mfrow=c(1,3), mar=c(2,1,2,1))
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$ref_time[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G1"]), col="darkgoldenrod1",
  stack=T, cex=.7, shrink=1, bins=100)
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$ref_time[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="S"]), col="coral1",
  stack=T, cex=.7, shrink=1, bins=100)
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$ref_time[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G2M"]), col="darkred",
  stack=T, cex=.7, shrink=1, bins=100)



# seurat prediction by seurat classes results
# plot(eval_seurat[[1]]$`5`$fit.seurat$ref_time,
#      eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift,
#      pch=16, cex=.7, xlab="")
library(circular)
par(mfrow=c(1,3), mar=c(2,1,2,1))
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G1"]), col="darkgoldenrod1",
  stack=T, cex=.7, shrink=1, bins=100)
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="S"]), col="coral1",
  stack=T, cex=.7, shrink=1, bins=100)
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G2M"]), col="darkred",
  stack=T, cex=.7, shrink=1, bins=100)




################################################

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



# compare cell time distribution
thetas <- lapply(1:length(inds), function(j) {
  ind <- inds[j]
  #  train_topX <- do.call(rbind, lapply(1:5, function(fold) {
  fl_name <- file.path("data/results/",
                       paste0("ind_",ind,"_data_withheld.rds"))
  df <- readRDS(fl_name)
  out <- df$theta.valid
  return(out)
})
names(thetas) <- inds

par(mfrow=c(2,3))
for (i in 1:length(thetas)) {
  hist(thetas[[i]], xlab="Fucci phase",
       main = inds[i], nclass=30)
}

# check in individual cyclical genes

fits_all <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
topgenes_all <- names(pve_all_ord)[1:5]


fits_NA18511 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA18511.rds")
fits_NA18855 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA18855.rds")
fits_NA18870 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA18870.rds")
fits_NA19098 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA19098.rds")
fits_NA19101 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA19101.rds")
fits_NA19160 <- readRDS("../output/trendfilter-individual.Rmd/fits_NA19160.rds")

pve_NA18511 <- sapply(fits_NA18511, "[[", 2)
pve_NA18855 <- sapply(fits_NA18855, "[[", 2)
pve_NA18870 <- sapply(fits_NA18870, "[[", 2)
pve_NA19098 <- sapply(fits_NA19098, "[[", 2)
pve_NA19101 <- sapply(fits_NA19101, "[[", 2)
pve_NA19160 <- sapply(fits_NA19160, "[[", 2)

pve_NA18511_ord <- pve_NA18511[order(pve_NA18511, decreasing = T)]
pve_NA18855_ord <- pve_NA18855[order(pve_NA18855, decreasing = T)]
pve_NA18870_ord <- pve_NA18870[order(pve_NA18870, decreasing = T)]
pve_NA19098_ord <- pve_NA19098[order(pve_NA19098, decreasing = T)]
pve_NA19101_ord <- pve_NA19101[order(pve_NA19101, decreasing = T)]
pve_NA19160_ord <- pve_NA19160[order(pve_NA19160, decreasing = T)]


topgenes <- list(names(pve_NA18511_ord)[1:5],
                 names(pve_NA18855_ord)[1:5],
                 names(pve_NA18870_ord)[1:5],
                 names(pve_NA19098_ord)[1:5],
                 names(pve_NA19101_ord)[1:5],
                 names(pve_NA19160_ord)[1:5])
names(topgenes) <- c("NA18511", "NA18855", "NA18870", "NA19098", "NA19101", "NA19160")

for (i in 1:length(topgenes)) {
  print(names(topgenes)[i])
  print(sum(topgenes_all %in% topgenes[[i]]))
}




