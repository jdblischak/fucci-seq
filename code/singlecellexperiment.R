#!/usr/bin/env Rscript

# Convert ExpressionSet to SingleCellExperiment class

library(Biobase)
library(SingleCellExperiment)

# Convert raw ------------------------------------------------------------------

eset_raw <- readRDS("data/eset-raw.rds")

# First convert to SummarizedExperiment
se_raw <- as(eset_raw, "SummarizedExperiment")

# Second convert to SingleCellExperiment
sce_raw <- as(se_raw, "SingleCellExperiment")

# Confirm identical
stopifnot(identical(exprs(eset_raw), assay(sce_raw)))
stopifnot(class(pData(eset_raw)) == "data.frame",
          class(colData(sce_raw)) == "DataFrame")
for (i in seq_len(ncol(pData(eset_raw)))) {
  stopifnot(identical(pData(eset_raw)[, i], colData(sce_raw)[, i]))
}

# Save
saveRDS(sce_raw, file = "data/sce-raw.rds")

# Convert final ----------------------------------------------------------------

eset_final <- readRDS("data/eset-final.rds")

# First convert to SummarizedExperiment
se_final <- as(eset_final, "SummarizedExperiment")

# Second convert to SingleCellExperiment
sce_final <- as(se_final, "SingleCellExperiment")

# Confirm identical
stopifnot(identical(exprs(eset_final), assay(sce_final)))
stopifnot(class(pData(eset_final)) == "data.frame",
          class(colData(sce_final)) == "DataFrame")
for (i in seq_len(ncol(pData(eset_final)))) {
  stopifnot(identical(pData(eset_final)[, i], colData(sce_final)[, i]))
}

# Save
saveRDS(sce_final, file = "data/sce-final.rds")
