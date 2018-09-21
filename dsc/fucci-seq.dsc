# project:
#   predicting cell cycle state using single-cell expression data
#
# description:
#   1. training prediction model, choose top X cyclical genes that optimizes
#      prediction error margin
#   2. compare results based on fucci cell times, fucci + dapi cell times,
#      predicting cells from mixed individuals, cells from one individual,
#   3. validation results applying the model to selecetd genes,
#      comparing the results with seurat method
#
# usage:
#


# Prepare data ------------------------------------------------------------------------

get_data_eset: R(eset <- readRDS(${fl_eset_raw}))
  $eset: eset

make_data_log2cpm: make_data_log2cpm.R + \
                   R(log2cpm <- make_data_log2cpm(eset))
  eset: $eset
  $log2cpm: log2cpm

make_data_normed: make_data_normed.R + \
                   R(log2cpm_normed <- make_data_normed(log2cpm, ncores))
  log2cpm: $logcpm
  $log2cpm_normed: log2cpm_normed
  ncores: 8

make_data_intensity: make_data_intensity.R + \
                   R(mat_intensity <- make_data_intensity(eset, type))
  eset: $eset
  type: fucci, fucci_dapi
  $mat_intensity: mat_intensity

make_data_theta: ${peco}/intensity2circle.R + \
                 R(theta <- intensity2circle(mat=mat_intensity, method) )
    mat_intensity: $mat_intensity
    method: trig, algebraic
    $theta: theta


# Predict cell times for cells from mixed individuals ----------------------------------

# prepare sample indices
make_data_partitions_individuals_mixed:
    make_data_partitions.R + \
    ${peco}\primes.R + ${peco}\partitionSamples.R + \
    R(partitions <- make_data_partitions(eset, training_valid=c(.85,.15), \
      training_folds=5, type="individual_mixed")) \
  eset: $eset
  $partitions: partitions

# model training: compute cyclical pve
training_get_cycle_pve.individuals_mixed: ${peco}/get.cyclical.R + \
                        R(which_training_sample <- partitions[[1]]$partitions[[fold]]$train
                          log2cpm_normed_train <- log2cpm_normed[,which_training_sample]
                          theta_train <- theta[,which_training_sample]
                          cyclical_pve <- get.cyclical(Y=log2cpm_normed_train,
                                                       theta=theta_train,
                                                       polyorder=2, ncores=ncores)
  log2cpm_normed: $log2com_normed
  fold: 1, 2, 3, 4, 5
  theta: $theta
  partitions: $partitions
  $cyclical_pve: cyclical_pve
  ncores: 15

# model training: prediction cell times in test samples
training_predict_test.individiuals_mixed:
  training_predict_test.R + \
  ${peco}/cycle.npreg.insample.R + \
  R(genes_topX <- rownames(cyclical_pve)[order(cyclical_pve,
                                          decreasing = T)[1:ngenes]];
    which_genes <- which(rownames(data_training$log2cpm.quant.nonvalid) %in% genes_topX)
    which_training_sample <- partitions[[1]]$partitions[[fold]]$train
    training_predict_test(log2cpm_normed, theta, pdata, which_training_sample,
                          which_test_sample, which_genes, ncores) )


# validation
validation.individuals_mixed:
  validation.R + \
  ${peco}/cycle.npreg.insample.R + \
  R(genes_topX <- rownames(cyclical_pve)[order(cyclical_pve,
                                          decreasing = T)[1:ngenes]];
    which_genes <- which(rownames(data_training$log2cpm.quant.nonvalid) %in% genes_topX)
    which_training_sample <- partitions[[1]]$partitions[[fold]]$train
    which_validation_sample <- partitions[[1]]$partitions[[fold]]$valid
    validation(log2cpm_normed, log2cpm,
                fdata, which_genes,
                which_training_sample,
                which_validation_sample,
                seurat.genes) )




# Predict cell times for cells from a particular individual ----------------------------


make_data_partitions_individuals_heldout: make_data_partitions.R + \
                      R(partitions <- make_data_partitions(eset, \
                                                           type="individual_heldout")) \
  eset: $eset
  $partitions: partitions

training_get_cycle_pve.individuals_heldout: ${peco}/get.cyclical.R + \
                        R(which_training_sample <- partitions[[ind]]$partitions[[fold]]$train
                          log2cpm_normed_train <- log2cpm_normed[,which_training_sample]
                          theta_train <- theta[,which_training_sample]
                          res <- get.cyclical(Y=log2cpm_normed_train,
                                              theta=theta_train,
                                              polyorder=2, ncores=ncores)
  log2cpm_normed: $log2com_normed
  fold: 1, 2, 3, 4, 5
  theta: $theta
  partitions: $partitions
  ncores: 15



#make_data_supervised:
#  make_data_supervised.R + \
#  R(pdata <- pData(eset); fdata <- fData(eset); \
#    data_list <- make_data_supervised(pdata, fdata, log2cpm, log2cpm_normed, theta))
#        eset: $eset
#        log2cpm: $log2cpm
#        log2cpm_normed: $log2com_normed
#        theta: $theta




DSC:
  R_libs: Biobase, parallel, stats, circular, conicfit
  define:
    make_data: get_data_eset * make_data_log2cpm * make_data_normed * make_data_intensity * make_data_theta * make_data_partitions * make_data_supervised

  run:
    first_pass: make_data

  exec_path: modules

  output: output

  global:
    dir_peco: "/project2/gilad/joycehsiao/fucci-seq/peco/R"
    fl_eset_raw: "/project2/gilad/joycehsiao/fucci-seq/data/eset-final.rds"

