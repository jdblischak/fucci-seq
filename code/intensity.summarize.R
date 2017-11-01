#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


library(EBImage)
library(methods)

# plates <- c("18511_18855", "18855_19101", "18855_19160", "18870_18511",
#             "18870_18855", "18870_19101", "18870_19160", "19098_18511",
#             "19098_18870", "19098_19160", "19101_18511", "19101_19098",
#             "19160_18870")
#args<-c("/project2/gilad/fucci-seq/intensities/","/project2/gilad/fucci-seq/intensities_stats/","18855_19101")
#args[1]=dir_intensity="/project2/gilad/fucci-seq/intensities/"
#args[2]=dir_intensity_stats="/project2/gilad/fucci-seq/intensities_stats/"
#args[3]=plate

#for (index in 1:length(plates)) {
#  plate <- plates[index]
  # extract foreground and background intensities of wells containing
  # a single nucleus
  nuclei <- readRDS(paste0(args[1], args[3], ".nuclei", ".rds"))

  singles <- names(nuclei)[which(unlist(nuclei) == 1)]

  # extract foreground/background intensities
  # of red and green channels

  #ints_fb <- do.call(rbind, lapply(1:5, function(index) {
  ints_fb <- do.call(rbind, lapply(1:length(singles), function(index) {
    id <- singles[index]
    ints <- readRDS(paste0(args[1], args[3], "/", args[3], ".", id, ".ints",".rds"))
    imgdata <- ints$imageOutput
    nnuclei <- nuclei[index]

    # ----- Approach 1: consider areas inside/outside the nucleus defined by DAPI
    # use EBImage functions to compute summary statistics

    # foreground intensity
    rfp.fore.stats <- computeFeatures.basic(imgdata$label,
                                imgdata$rfp, basic.quantiles = c(.1, .5, .9))
    gfp.fore.stats <- computeFeatures.basic(imgdata$label,
                                 imgdata$gfp, basic.quantiles = c(.1, .5, .9))
    dapi.fore.stats <- computeFeatures.basic(imgdata$label,
                                 imgdata$dapi, basic.quantiles = c(.1, .5, .9))
    size <- computeFeatures.shape(imgdata$label, imgdata$dapi)[,1]

    #background intensity
    rfp.back.stats <- computeFeatures.basic(1-imgdata$label,
                                       imgdata$rfp, basic.quantiles = c(.1, .5, .9))
    gfp.back.stats <- computeFeatures.basic(1-imgdata$label,
                                       imgdata$gfp, basic.quantiles = c(.1, .5, .9))
    dapi.back.stats <- computeFeatures.basic(1-imgdata$label,
                                        imgdata$dapi, basic.quantiles = c(.1, .5, .9))

    # rename *.stats
    dimnames(rfp.fore.stats)[[2]] <- gsub("b.", "rfp.fore.", dimnames(rfp.fore.stats)[[2]])
    dimnames(gfp.fore.stats)[[2]] <- gsub("b.", "gfp.fore.", dimnames(gfp.fore.stats)[[2]])
    dimnames(dapi.fore.stats)[[2]] <- gsub("b.", "dapi.fore.", dimnames(dapi.fore.stats)[[2]])
    dimnames(rfp.back.stats)[[2]] <- gsub("b.", "rfp.back.", dimnames(rfp.back.stats)[[2]])
    dimnames(gfp.back.stats)[[2]] <- gsub("b.", "gfp.back.", dimnames(gfp.back.stats)[[2]])
    dimnames(dapi.back.stats)[[2]] <- gsub("b.", "dapi.back.", dimnames(dapi.back.stats)[[2]])


    # ----- Approach 2: consider the 100 x 100 pixel intensity area centered at the nucleus
    # note that the area is bounded by the distance of the center to the image edges

    # foreground intensity
    rfp.fore.stats.zoom <- computeFeatures.basic(imgdata$label.zoom, imgdata$rfp)
    gfp.fore.stats.zoom <- computeFeatures.basic(imgdata$label.zoom, imgdata$gfp)
    dapi.fore.stats.zoom <- computeFeatures.basic(imgdata$label.zoom, imgdata$dapi)
    size.zoom <- computeFeatures.shape(imgdata$label.zoom, imgdata$dapi)[,1]

    #background intensity
    rfp.back.stats.zoom <- computeFeatures.basic(1-imgdata$label.zoom, imgdata$rfp)
    gfp.back.stats.zoom <- computeFeatures.basic(1-imgdata$label.zoom, imgdata$gfp)
    dapi.back.stats.zoom <- computeFeatures.basic(1-imgdata$label.zoom, imgdata$dapi)

    # summary stats
    rfp.sum.zoom <- (rfp.fore.stats.zoom$b.mean-rfp.back.stats.zoom$b.mean)*size.zoom
    gfp.sum.zoom <- (gfp.fore.stats.zoom$b.mean-gfp.back.stats.zoom$b.mean)*size.zoom
    dapi.sum.zoom <- (dapi.fore.stats.zoom$b.mean-dapi.back.stats.zoom$b.mean)*size.zoom

    # rename *.stats
    dimnames(rfp.fore.stats)[[2]] <- gsub("b.", "rfp.fore.", dimnames(rfp.fore.stats)[[2]])
    dimnames(gfp.fore.stats)[[2]] <- gsub("b.", "gfp.fore.", dimnames(gfp.fore.stats)[[2]])
    dimnames(dapi.fore.stats)[[2]] <- gsub("b.", "dapi.fore.", dimnames(dapi.fore.stats)[[2]])
    dimnames(rfp.back.stats)[[2]] <- gsub("b.", "rfp.back.", dimnames(rfp.back.stats)[[2]])
    dimnames(gfp.back.stats)[[2]] <- gsub("b.", "gfp.back.", dimnames(gfp.back.stats)[[2]])
    dimnames(dapi.back.stats)[[2]] <- gsub("b.", "dapi.back.", dimnames(dapi.back.stats)[[2]])


    data.frame(wellID=id, size=size, nnuclei=nnuclei,
               rfp.fore.stats, gfp.fore.stats, dapi.fore.stats,
               rfp.back.stats, gfp.back.stats, dapi.back.stats,
               rfp.sum.zoom=rfp.sum.zoom, gfp.sum.zoom=gfp.sum.zoom, dapi.sum.zoom=dapi.sum.zoom)
  }) )

  saveRDS(ints_fb, file = paste0(args[2], args[3], ".stats.rds"))
