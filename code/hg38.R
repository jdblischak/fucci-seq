#!/usr/bin/env Rscript

# Hacky script to compare hg19 molecule counts in ~/fucci-seq/data/molecules and
# hg38 molecules counts in ~/fucci-seq-hg38/data/molecules.
#
# Script executed in ~/fucci-seq-hg38/

library(data.table)
library(stringr)

# Input ------------------------------------------------------------------------

mol_hg19 <- fread("zcat ../fucci-seq/data/molecules/*.txt.gz | grep -v sample")
colnames(mol_hg19) <- colnames(fread("zcat ../fucci-seq/data/molecules/20170905.txt.gz",
                                     nrows = 1))
dim(mol_hg19)
mol_hg19[, 1:5]

mol_hg38 <- fread("zcat data/molecules/*.txt.gz | grep -v sample")
colnames(mol_hg38) <- colnames(fread("zcat data/molecules/20170905.txt.gz",
                                     nrows = 1))
dim(mol_hg38)
mol_hg38[, 1:5]

# The sample should be identical
stopifnot(mol_hg19$sample == mol_hg38$sample,
          mol_hg19$experiment == mol_hg38$experiment,
          mol_hg19$well == mol_hg38$well)

# There are slight differences in annotation between releases 75 and 96
genes_hg19 <- str_subset(colnames(mol_hg19), "^ENSG")
genes_hg38 <- str_subset(colnames(mol_hg38), "^ENSG")
genes_all <- union(genes_hg19, genes_hg38)
genes_shared <- intersect(genes_hg19, genes_hg38)
length(genes_shared) / length(genes_all)
## [1] 0.8929394

# Limiting to protein-coding genes
saf <- fread("/project2/gilad/fucci-seq/genome-ensembl-release-75/hs.saf")
saf
protein_coding <- unique(saf$GeneID)
genes_hg19 <- genes_hg19[genes_hg19 %in% protein_coding]
genes_hg38 <- genes_hg38[genes_hg38 %in% protein_coding]
genes_all <- union(genes_hg19, genes_hg38)
genes_shared <- intersect(genes_hg19, genes_hg38)
length(genes_shared) / length(genes_all)
## [1] 0.9351109

# Filter -----------------------------------------------------------------------

cols_shared <- c("sample", "experiment", "well", genes_shared)

mol_hg19_shared <- mol_hg19[, ..cols_shared]
dim(mol_hg19_shared)
mol_hg38_shared <- mol_hg38[, ..cols_shared]
dim(mol_hg38_shared)
stopifnot(dim(mol_hg19_shared) == dim(mol_hg38_shared),
          colnames(mol_hg19_shared) == colnames(mol_hg38_shared))


# Remove zeros
sums_hg19 <- mol_hg19_shared[, vapply(.SD, sum, integer(1)), .SDcols = -(sample:well)]
sums_hg38 <- mol_hg38_shared[, vapply(.SD, sum, integer(1)), .SDcols = -(sample:well)]
cor(sums_hg19, sums_hg38)
plot(sums_hg19, sums_hg38)
sums_total <- sums_hg19 + sums_hg38
zeros <- sums_total == 0
sum(zeros)
genes_shared_nonzero <- genes_shared[!zeros]
cols_shared_nonzero <- c("sample", "experiment", "well", genes_shared_nonzero)
mol_hg19_shared_nonzero <- mol_hg19_shared[, ..cols_shared_nonzero]
dim(mol_hg19_shared_nonzero)
mol_hg38_shared_nonzero <- mol_hg38_shared[, ..cols_shared_nonzero]
dim(mol_hg38_shared_nonzero)
stopifnot(dim(mol_hg19_shared_nonzero) == dim(mol_hg38_shared_nonzero),
          colnames(mol_hg19_shared_nonzero) == colnames(mol_hg38_shared_nonzero))

# Correlation ------------------------------------------------------------------

corrs <- numeric(length = length(genes_shared_nonzero))
names(corrs) <- genes_shared_nonzero
for (i in seq_along(genes_shared_nonzero)) {
  counts_hg19 <- unlist(mol_hg19_shared_nonzero[, genes_shared_nonzero[i], with = FALSE])
  counts_hg38 <- unlist(mol_hg38_shared_nonzero[, genes_shared_nonzero[i], with = FALSE])
  corrs[i] <- cor(counts_hg19, counts_hg38)
}
summary(corrs)
hist(corrs)

plot(unlist(mol_hg19_shared_nonzero[, 5]), unlist(mol_hg38_shared_nonzero[, 5]))
cor(unlist(mol_hg19_shared_nonzero[, 5]), unlist(mol_hg38_shared_nonzero[, 5]))

# Plot -------------------------------------------------------------------------
