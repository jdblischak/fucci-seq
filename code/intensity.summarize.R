# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
#

library(EBImage)
library(methods)

plates <- c("18511_18855", "18855_19101", "18855_19160", "18870_18511",
            "18870_18855", "18870_19101", "18870_19160", "19098_18511",
            "19098_18870", "19098_19160", "19101_18511", "19101_19098",
            "19160_18870")

dir_intensity <- "/project2/gilad/fucci-seq/intensities/"
dir_intensity_stats  <- "/project2/gilad/fucci-seq/intensities_stats/"

for (index in 1:length(plates)) {
  plate <- plates[index]
  # extract foreground and background intensities of wells containing
  # a single nucleus
  nuclei <- readRDS(paste0(dir_intensity, plate, ".nuclei", ".rds"))
  singles <- names(nuclei)[which(unlist(nuclei) == 1)]

  # extract foreground/background intensities
  # of red and green channels

  #ints_fb <- do.call(rbind, lapply(1:5, function(index) {
  ints_fb <- do.call(rbind, lapply(1:length(singles), function(index) {
    id <- singles[index]
    ints <- readRDS(paste0(dir_intensity, plate, "/", plate, ".", id, ".ints",".rds"))
    imgdata <- ints$imageOutput

    data.frame(rfp.fore.mean = mean(imgdata$rfp[imgdata$label == 1]),
               rfp.back.mean = mean(imgdata$rfp[imgdata$label == 0]),
               gfp.fore.mean = mean(imgdata$gfp[imgdata$label == 1]),
               gfp.back.mean = mean(imgdata$gfp[imgdata$label == 0]),
              rfp.fore.z.mean = mean(scale(imgdata$rfp)[imgdata$label == 1]),
              rfp.back.z.mean = mean(scale(imgdata$rfp)[imgdata$label == 0]),
              gfp.fore.z.mean = mean(scale(imgdata$gfp)[imgdata$label == 1]),
              gfp.back.z.mean = mean(scale(imgdata$gfp)[imgdata$label == 0]))
  }) )
  ints_fb$well <- singles

  saveRDS(ints_fb, file = paste0(dir_intensity_stats, plate, ".stats.rds"))
}
