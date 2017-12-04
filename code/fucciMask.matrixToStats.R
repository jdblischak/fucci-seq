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
    dapi.fore <- obj$dapi[obj$label==1]
    dapi.back <- obj$dapi[obj$label==0]

    rfp.backmean.log.sum <- log10(sum(rfp.fore-mean(rfp.back)))
    rfp.backmedian.log.sum <- log10(sum(rfp.fore-median(rfp.back)))

    gfp.backmean.log.sum <- log10(sum(gfp.fore-mean(gfp.back)))
    gfp.backmedian.log.sum <- log10(sum(gfp.fore-median(gfp.back)))

    dapi.backmean.log.sum <- log10(sum(dapi.fore-mean(dapi.back)))
    dapi.backmedian.log.sum <- log10(sum(dapi.fore-median(dapi.back)))

    rm(obj)
    data.frame(plate=plate, id=id,
               rfp.backmean.log.sum, gfp.backmean.log.sum, dapi.backmean.log.sum,
               rfp.backmedian.log.sum, gfp.backmedian.log.sum, dapi.backmedian.log.sum)
  }))
  return(res)
}))

saveRDS(ints, file="/project2/gilad/fucci-seq/intensities_matrix2stats/ints.rds")



