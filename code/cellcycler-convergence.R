#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Description:
#   This file runs parametric cellcycleR at a fixed random seed.
#
# Usage:
#
# Rscript cellcycler-convergence.R ${seed} ${iter}


args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}



# ----- function starts here

library("cellcycleR")

ints <- readRDS(file="/project2/gilad/joycehsiao/fucci-seq/data/intensity.rds")

ints$rfp.z <- scale(ints$rfp.median.log10sum, center=TRUE, scale=TRUE)
ints$gfp.z <- scale(ints$gfp.median.log10sum, center=TRUE, scale=TRUE)
ints$dapi.z <- scale(ints$dapi.median.log10sum, center=TRUE, scale=TRUE)

ints_tmp <- data.frame(rfp.z=ints$rfp.z,
                       gfp.z=ints$gfp.z,
                       dapi.z=ints$dapi.z)

set.seed(args[1])
sessionInfo()

out <- sin_cell_ordering_class(as.matrix(ints_tmp),
                               celltime_levels=200,
                               num_shuffle=1, tol=1e-6,
                               maxiter=500, n_cores=8)

saveRDS(out,
        file = paste0("/project2/gilad/joycehsiao/fucci-seq/output_tmp/cellcycer-convergence/out-",
                      args[2],".rds"))

