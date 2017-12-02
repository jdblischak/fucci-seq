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

#' @example
# plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
#             "18870_18855","18870_19101","18870_19160","19098_18511",
#             "19098_18870","19098_19160","19101_18511","19101_19098",
#             "19160_18870","19101_19160","19160_18511", "18855_19098")
# plates <- c("18855_19101","18855_19160","18870_18511",
#             "18870_18855","18870_19101","18870_19160","19098_18511",
#             "19098_18870","19098_19160","19101_18511","19101_19098",
#             "19160_18870","19101_19160","19160_18511", "18855_19098")
# for (index in 1:length(plates)) {
#   fucciMask.convertToMatrix(plates[index])
# }

#plate <- args[1]

fucciMask.convertToMatrix <- function (plate) {

  library(EBImage)
  message("converting images for ", plate)

  dir_ints <- paste0("/project2/gilad/fucci-seq/intensities")

  # output directory
  dir_output <- "/project2/gilad/joycehsiao/fucci-seq/output_tmp"

  # read in nuclei information for the plate
  fl_nuclei <- paste0(plate, ".nuclei.rds")
  df_nuclei <- readRDS(file = file.path(dir_ints, fl_nuclei))

  # read in intensity data
  fl_list <- list.files(file.path(dir_ints, plate), full.names = TRUE)
  fl_well_labels <- sapply(strsplit(basename(fl_list), split=".", fixed = TRUE), "[[", 2)

  # select intensity files that have one nucleus
  singles <- names(df_nuclei)[which(df_nuclei==1)]
  assertthat::are_equal(fl_well_labels[which(fl_well_labels %in% singles)], singles)

  fl_singles <- fl_list[which(fl_well_labels %in% singles)]

  # for every well, convert EBImage object to matrix
  for (index in 1:length(fl_singles)) {
      message("well ", singles[index], "\n")

      df_ints <- readRDS(file = fl_singles[index])

      # convert image data to matrices
      # code this into the part where we save the image intensity files
      label <- df_ints$imageOutput$label
      label <- label@.Data

      label.zoom <- df_ints$imageOutput$label.zoom
      label.zoom <- label.zoom@.Data

      rfp <- df_ints$imageOutput$rfp
      rfp <- rfp@.Data

      gfp <- df_ints$imageOutput$gfp
      gfp <- gfp@.Data

      ints_tmp <- data.frame(label=c(label), label.zoom=c(label.zoom),
                         rfp = c(rfp), gfp = c(gfp), dapi=dapi())
      saveRDS(ints_tmp,
              file = paste0(dir_output, "/", plate, ".", singles[index],".tmp.rds"))
  }

  message("combine well data to one rds", "\n")

  # read RDS
  ints_long <- do.call(rbind, lapply(1:length(singles), function(index) {
    ints <- readRDS(paste0(dir_output, "/", plate, ".", singles[index],".tmp.rds"))
    ints$well <- singles[index]
    return(ints)
  }) )
  saveRDS(ints_long, paste0(dir_output, "/", plate, ".tmp.long.rds"))

}

fucciMask.convertToMatrix(args[1])


