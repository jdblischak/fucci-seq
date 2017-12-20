#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Description:
#   This file runs parametric cellcycleR at a fixed random seed.
#
# Usage:
#
# Rscript mclapply.R ${seed} ${iter}


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

#n_cores <- parallel::detectCores()
#n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

set.seed(args[1])
# use lapply instead of mclapply
#out <- sin_cell_ordering_class(as.matrix(ints_tmp),
#                               celltime_levels=200, num_iter=10,
#			       verbose=FALSE,freq=1,
#			       fix.phase=FALSE, phase_in=NULL,
#			       parallel = FALSE,  n_cores=NULL)

# use mclapply
out <- system.time(sin_cell_ordering_class(as.matrix(ints_tmp),
                               celltime_levels=200, num_iter=100,
                               fix.phase=FALSE, phase_in=NULL,
			       freq=1, verbose=FALSE,
                               parallel = TRUE,  n_cores=6))

saveRDS(out,
        file = paste0("/project2/gilad/joycehsiao/fucci-seq/code/scratch/mclapply-G3.rds"))

