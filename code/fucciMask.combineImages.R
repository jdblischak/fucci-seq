#!/usr/bin/env Rscript

# Command line arguments can be read with commandArgs.
#
# Description:
#   This file contains an example of how to run through the image
# processing code. It runs through all 96 images, saves the data frame,
# then builds up an image containing zoomed in, aligned images of all
# identified nuclei in all channels.
#
# Usage:
# dir_code="/project2/gilad/joycehsiao/fucci-seq/code/"
# dir_images_wellId="/project2/gilad/fucci-seq/images_wellid/"
# dir_intensities="/project2/gilad/fucci-seq/intensities/"
# dir_images_processed="/scratch/midway2/joycehsiao/fucci-seq/images_processed/18511_18855/"
# dir_images_combined="/project2/gilad/fucci-seq/images_combined/"
# pl="18511_18855"
#
# Rscript ${dir_code}/fucciMask.combineImages.R ${dir_images_wellId} ${dir_intensities} ${dir_images_processed} ${dir_images_combined} ${pl}
#

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}

# specify arguments
#dir_images_wellId <- args[1]
#dir_intensities <- args[2]
#dir_images_processed <- args[3]
#dir_images_combined <- args[4]
#pl <- args[5]


# ----- function starts here

library("EBImage")
library("methods")

# import results on the number of nuclei
nuclei <- unlist(readRDS(paste0(args[2], args[5], ".nuclei", ".rds")))

# grab all the rows with exactly one nucleus
singles <- names(nuclei)[which(nuclei == 1)]

# make a nice figure. This code is a terrible mess. It creates red,
# green, and blue images for every row/column. I show bg subtracted
# data for red and green channels.

red.0 <- list()
green.0 <- list()
blue.0 <- list()

red.1   <- list()
green.1 <- list()
blue.1  <- list()

red.2   <- list()
green.2 <- list()
blue.2  <- list()

red.3   <- list()
green.3 <- list()
blue.3  <- list()

red.4   <- list()
green.4 <- list()
blue.4  <- list()

for (i in 1:length(singles))
{
  id <- singles[i]

  # ---- prepare for well number
  well <- readImage(paste0( args[1], id, ".tiff"))
  red.0[[i]]   <- well[,,1]
  green.0[[i]] <- well[,,2]
  blue.0[[i]]  <- well[,,3]

  # --- prepare for channel images
  bright <- readImage(paste0(args[3], id, '.1.bright.TIFF'))
  nucs   <- readImage(paste0(args[3], id, '.1.outline.TIFF'))
  gfps   <- readImage(paste0(args[3], id, '.1.gfp.TIFF'))
  rfps   <- readImage(paste0(args[3], id, '.1.rfp.TIFF'))

  # convert from images to 3D matrices
  bright <- imageData(rgbImage(bright, bright, bright))
  # to make greyscale, red=blue=green
  nucs   <- imageData(nucs)
  gfps   <- imageData(rgbImage(green=gfps))
  rfps   <- imageData(rgbImage(red=rfps))

  red.1[[i]]   <- bright[,,1]
  green.1[[i]] <- bright[,,2]
  blue.1[[i]]  <- bright[,,3]

  red.2[[i]]   <- nucs[,,1]
  green.2[[i]] <- nucs[,,2]
  blue.2[[i]]  <- nucs[,,3]

  red.3[[i]]   <- gfps[,,1]
#  green.3[[i]] <- gfps[,,2] - singles$bg.gfp[i] # subtract background
  green.3[[i]] <- gfps[,,2]
  blue.3[[i]]  <- gfps[,,3]

#  red.4[[i]]   <- rfps[,,1] - singles$bg.rfp[i] # subtract background
  red.4[[i]]   <- rfps[,,1]
  green.4[[i]] <- rfps[,,2]
  blue.4[[i]]  <- rfps[,,3]
}

reds.0   <- as.Image(do.call('cbind', red.0))
greens.0 <- as.Image(do.call('cbind', green.0))
blues.0  <- as.Image(do.call('cbind', blue.0))

reds.1   <- as.Image(do.call('cbind', red.1))
greens.1 <- as.Image(do.call('cbind', green.1))
blues.1  <- as.Image(do.call('cbind', blue.1))

reds.2   <- as.Image(do.call('cbind', red.2))
greens.2 <- as.Image(do.call('cbind', green.2))
blues.2  <- as.Image(do.call('cbind', blue.2))

reds.3   <- as.Image(do.call('cbind', red.3))
greens.3 <- normalize(as.Image(do.call('cbind', green.3)))
# normalize to bring scale back up to 0-1 after bg removal
blues.3  <- as.Image(do.call('cbind', blue.3))

reds.4   <- normalize(as.Image(do.call('cbind', red.4)))
# normalize to bring scale back up to 0-1 after bg removal
greens.4 <- as.Image(do.call('cbind', green.4))
blues.4  <- as.Image(do.call('cbind', blue.4))

# combine images from well ID, and the four channels
combined <- rgbImage(
  as.Image(rbind(imageData(reds.0),
          imageData(reds.1),imageData(reds.2),imageData(reds.3),imageData(reds.4))),
  as.Image(rbind(imageData(greens.0),
          imageData(greens.1),imageData(greens.2),imageData(greens.3),imageData(greens.4))),
  as.Image(rbind(imageData(blues.0),
          imageData(blues.1),imageData(blues.2),imageData(blues.3),imageData(blues.4))))

writeImage(combined, file=paste0(args[4], args[5], ".combined.TIFF"))
