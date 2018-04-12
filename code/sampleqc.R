#!/usr/bin/env Rscript

# Utility script to combine the ExpressionSet RDS files for each chip into a
# combined ExpressionSet RDS object for quick import.
#
# Usage:
#
# Rscript output-combined-eset.R <directory with ExpressionSet rds files> \
#                                <output combined RDS file>
#

suppressPackageStartupMessages(library("Biobase"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dir_eset <- args[1]
stopifnot(dir.exists(dir_eset))
outfile <- args[2]

fname <- Sys.glob(file.path(dir_eset, "*.rds"))
eset <- Reduce(combine, Map(readRDS, fname))

saveRDS(eset, file = outfile)




## Prepare data
library("Biobase")
fname <- Sys.glob("../data/eset/*.rds")
eset <- Reduce(combine, Map(readRDS, fname))
anno <- pData(eset)

## Total mapped reads cutoff
cut_off_reads <- quantile(anno[anno$cell_number == 0,"mapped"], 0.82)
anno$cut_off_reads <- anno$mapped > cut_off_reads

## Unmapped ratio cutoff
anno$unmapped_ratios <- anno$unmapped/anno$umi
cut_off_unmapped <- quantile(anno[anno$cell_number == 0,"unmapped_ratios"], 0.40)
anno$cut_off_unmapped <- anno$unmapped_ratios < cut_off_unmapped

## ERCC percentage cutoff
anno$ercc_percentage <- anno$reads_ercc / anno$mapped
cut_off_ercc <- quantile(anno[anno$cell_number == 0,"ercc_percentage"], 0.20)
anno$cut_off_ercc <- anno$ercc_percentage < cut_off_ercc

## Number of genes detected cutoff
cut_off_genes <- quantile(anno[anno$cell_number == 0,"detect_hs"], 0.80)
anno$cut_off_genes <- anno$detect_hs > cut_off_genes

## Total molecule outlier
## create 3 groups according to cell number
group_3 <- rep("two",dim(anno)[1])
group_3[grep("0", anno$cell_number)] <- "no"
group_3[grep("1", anno$cell_number)] <- "one"
## create data frame
library(dplyr)
data <- anno %>% dplyr::select(experiment:concentration, mapped, molecules)
data <- data.frame(data, group = group_3)

## perform lda
library("MASS")
data_lda <- lda(group ~ concentration + molecules, data = data)
data_lda_p <- predict(data_lda, newdata = data[,c("concentration", "molecules")])$class
data$data_lda_p <- data_lda_p
## identify the outlier
library("tibble")
outliers_lda <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_lda_p == "two")
## create filter
anno$molecule_outlier <- row.names(anno) %in% outliers_lda$sample_id

## Read to molecule conversion outlier
## calculate convertion
anno$ercc_conversion <- anno$mol_ercc / anno$reads_ercc
anno$conversion <- anno$mol_hs / anno$reads_hs
## try lda
data$conversion <- anno$conversion
data$ercc_conversion <- anno$ercc_conversion
data_ercc_lda <- lda(group ~ ercc_conversion + conversion, data = data)
data_ercc_lda_p <- predict(data_ercc_lda,  newdata = data[,c("ercc_conversion", "conversion")])$class
data$data_ercc_lda_p <- data_ercc_lda_p
## identify the outlier
outliers_conversion <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_ercc_lda_p == "two")
## create filter
anno$conversion_outlier <- row.names(anno) %in% outliers_conversion$sample_id

## Final filter
anno$filter_all <- anno$cell_number == 1 &
                   anno$mol_egfp > 0 &
                   anno$valid_id &
                   anno$cut_off_reads &
                   anno$cut_off_unmapped &
                   anno$cut_off_ercc &
                   anno$cut_off_genes &
                   anno$molecule_outlier == "FALSE" &
                   anno$conversion_outlier == "FALSE"
