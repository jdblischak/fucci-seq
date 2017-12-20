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

G <- 500;
num_cells <- 400
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);

celltime_levels <- 100;

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
out <- system.time(sin_cell_ordering_class(cycle_data,
                               celltime_levels = 100, num_iter=500,
                               fix.phase=FALSE, phase_in=NULL,
                               freq=1, verbose=FALSE,
                               parallel=TRUE, n_cores=8))

saveRDS(out,
        file = paste0("/project2/gilad/joycehsiao/fucci-seq/code/scratch/mclapply-G100.rds"))

