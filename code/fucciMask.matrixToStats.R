#!/usr/bin/env Rscript

# Description:
# This script takes the output rds from fucciMas.R
# and converts the pixel intensity matrix from EBImage object to matrix
#
# The main reason is that the current version of EBImage cannot be
# installed directly on broadwl or mstephens partition, so I can't
# use them in rstudio cloud...
#
# To run this script on the command line

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}



 plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
             "18870_18855","18870_19101","18870_19160","19098_18511",
             "19098_18870","19098_19160","19101_18511","19101_19098",
             "19160_18870","19101_19160","19160_18511", "18855_19098")

#plates <- c("18511_18855","18855_19101")

ints <- do.call(rbind, lapply(1:length(plates), function(index) {
  plate <- plates[index]
  nuclei <- readRDS(paste0("/project2/gilad/fucci-seq/intensities/",plate,".nuclei.rds"))
  singles <- nuclei[nuclei==1]
  singles <- names(singles)[!is.na(singles)]

  res <- do.call(rbind, lapply(1:length(singles), function(index2){
    id <- singles[index2]
    obj <- readRDS(paste0("/project2/gilad/joycehsiao/fucci-seq/output_tmp/",
                          plate,".",id,".tmp.rds"))
    rfp.fore <- obj$rfp[obj$label.zoom==1]
    rfp.back <- obj$rfp[obj$label.zoom==0]
    gfp.fore <- obj$gfp[obj$label.zoom==1]
    gfp.back <- obj$gfp[obj$label.zoom==0]
    rfp.backmean.sum.log <- log(sum(rfp.fore-mean(rfp.back)))
    gfp.backmean.sum.log <- log(sum(gfp.fore-mean(gfp.back)))
    rfp.backmedian.sum.log <- log(sum(rfp.fore-median(rfp.back)))
    gfp.backmedian.sum.log <- log(sum(gfp.fore-median(gfp.back)))
    rm(obj)
    data.frame(plate=plate, id=id, rfp.backmean.sum.log, gfp.backmean.sum.log,
               rfp.backmedian.sum.log, gfp.backmedian.sum.log)
  }))
  return(res)
}))

saveRDS(ints, file="/project2/gilad/fucci-seq/intensities_matrix2stats/ints.rds")



