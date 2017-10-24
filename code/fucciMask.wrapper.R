#!/usr/bin/env Rscript

# Usage:
# dir_images_data="/project2/gilad/fucci-seq/images_curated/18511_18855/"
# dir_code="/project2/gilad/joycehsiao/fucci-seq/code/"
# dir_images_processed="/scratch/midway2/joycehsiao/fucci-seq/images_processed/18511_18855/"
# dir_intensities="/project2/gilad/fucci-seq/intensities/"
# pl="18511_18855"
#
# args <- c(dir_images_data, dir_code, dir_images_processed, dir_intensities, pl)
# Rscript ./fucciMask.wrapper.R ${dir_images_data} ${dir_code} ${dir_images_processed} ${dir_intensities} ${pl}
#
#' @title fucciMask wrapper
#' @description This code runs through all 96 images, saves the data frame,
# then builds up an image containing zoomed in, aligned images of all
# identified nuclei in all channels.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

library('EBImage')
source(paste0(args[2], "fucciMask.R"))

# run through and read in all the images
data <- vector("list", 96)
for (i in 1:96)
{
  id <- sprintf("%05d", i)
  message("Processing image ", id, "\n")
  message("Reading images...")
  bright <- readImage(paste0(args[1], "BRIGHT/", id, ".TIFF"))
  dapi   <- readImage(paste0(args[1], "DAPI/",  id, ".TIFF"))
  gfp    <- readImage(paste0(args[1], "GFP/",   id,'.TIFF'))
  rfp    <- readImage(paste0(args[1], "RFP/",   id,'.TIFF'))

  message("Analyzing image...")
  tmp <- create_mask(bright, dapi, gfp, rfp, id, dir_output=args[3],
                      control=list(medianFilterRadius=10,
                                   size_cutoff=350, display=FALSE,
                                   printWholeDAPI=FALSE,
                                   printProcessedImages=TRUE,
                                   displayProcessedImages=FALSE))
  data[[i]] <- tmp
  saveRDS(data[[i]],
    file = paste0( args[4], args[5],"/", args[5],".", id, ".ints.rds"))
}
#names(data) <- sapply(data, "[[", 1)
#saveRDS(data, file = paste0(args[4], "ints_", args[5] ,".rds"))

# save nuclei count in a separate rds
nuclei_num <- lapply(data, "[[", 2)
names(nuclei_num) <- sapply(data, "[[", 1)
nuclei_num <- do.call(c, nuclei_num)
saveRDS(nuclei_num, file = paste0(args[4], args[5] ,".nuclei",".rds"))
