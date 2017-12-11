# Combine intensity data from different plates
#
# Output: fucci-seq/data/intensity.rds
#
# Description:
#   Combine intensity summary data from different plates, and
# also add individual ID data from fucci-seq/data/eset. Output
# is used for exploratory data analysis purposes. Later, the
# intensity data will be combined with the eset data.
#
# Date updated: 2017-12-11

plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
            "18870_18855","18870_19101","18870_19160","19098_18511",
            "19098_18870","19098_19160","19101_18511","19101_19098",
            "19160_18870","19101_19160","19160_18511", "18855_19098")

##Combine intensity stats from different plates.
ints <- do.call(rbind, lapply(1:length(plates), function(index) {
  tmp <- readRDS(paste0("/project2/gilad/fucci-seq/intensities_stats/",plates[index],".stats.rds"))
  tmp <- data.frame(plate=plates[index],
                    image=as.character(droplevels(tmp$wellID)),
                    size=tmp$size,
                    perimeter=tmp$perimeter,
                    eccentricity=tmp$eccentricity,
                    rfp.fore.zoom.mean=tmp$rfp.fore.zoom.mean,
                    rfp.fore.zoom.median=tmp$rfp.fore.zoom.q05,
                    gfp.fore.zoom.mean=tmp$gfp.fore.zoom.mean,
                    gfp.fore.zoom.median=tmp$gfp.fore.zoom.q05,
                    dapi.fore.zoom.mean=tmp$dapi.fore.mean,
                    dapi.fore.zoom.median=tmp$dapi.fore.q05,
                    rfp.back.zoom.mean=tmp$rfp.back.zoom.mean,
                    rfp.back.zoom.median=tmp$rfp.back.zoom.q05,
                    gfp.back.zoom.mean=tmp$gfp.back.zoom.mean,
                    gfp.back.zoom.median=tmp$gfp.back.zoom.q05,
                    dapi.back.zoom.mean=tmp$dapi.back.mean,
                    dapi.back.zoom.median=tmp$dapi.back.q05,
                    rfp.mean.log10sum=with(tmp, log10((rfp.fore.zoom.mean-rfp.back.zoom.mean)*size.zoom)),
                    gfp.mean.log10sum=with(tmp, log10((gfp.fore.zoom.mean-gfp.back.zoom.mean)*size.zoom)),
                    dapi.mean.log10sum=with(tmp, log10((dapi.fore.mean-dapi.back.mean)*size)),
                    rfp.median.log10sum=with(tmp, log10((rfp.fore.zoom.mean-rfp.back.zoom.q05)*size.zoom)),
                    gfp.median.log10sum=with(tmp, log10((gfp.fore.zoom.mean-gfp.back.zoom.q05)*size.zoom)),
                    dapi.median.log10sum=with(tmp, log10((dapi.fore.mean-dapi.back.q05)*size)) )

  tmp$rfp.mean.log10sum[which(tmp$rfp.mean.log10sum == "NaN")] <- min(tmp$rfp.mean.log10sum, na.rm=TRUE)
  tmp$gfp.mean.log10sum[which(tmp$gfp.mean.log10sum == "NaN")] <- min(tmp$gfp.mean.log10sum, na.rm=TRUE)
  tmp$dapi.mean.log10sum[which(tmp$dapi.mean.log10sum == "NaN")] <- min(tmp$dapi.mean.log10sum, na.rm=TRUE)
  tmp$rfp.median.log10sum[which(tmp$rfp.median.log10sum == "NaN")] <- min(tmp$rfp.median.log10sum, na.rm=TRUE)
  tmp$gfp.median.log10sum[which(tmp$gfp.median.log10sum == "NaN")] <- min(tmp$gfp.median.log10sum, na.rm=TRUE)
  tmp$dapi.median.log10sum[which(tmp$dapi.median.log10sum == "NaN")] <- min(tmp$dapi.median.log10sum, na.rm=TRUE)

  return(tmp)
}) )


##Load the expression set info.
eset_fls <- list.files("/project2/gilad/joycehsiao/fucci-seq/data/eset", full.names=TRUE)

anno <- do.call(rbind, lapply(1:length(eset_fls), function(index) {
  eset_index <- readRDS(eset_fls[[index]])
  pdata_index <- pData(eset_index)
  return(pdata_index)
}))

# make unique id in both
ints$unique <- paste0(ints$plate,"_",as.character(ints$image))
anno$unique <- paste0(anno$image_individual,"_", sprintf("%05d",as.numeric(anno$image_label)))

# subset annotation file
subset_index <- which(anno$unique %in% ints$unique)
anno_subset <- anno[subset_index,]

# subset intensities
ints_subset_index <- match(anno_subset$unique, ints$unique)
ints_subset <- ints[ints_subset_index,]

#all.equal(anno_subset$unique, ints_subset$unique)
ints_subset$chip_id <- anno_subset$chip_id


saveRDS(ints_subset, file="/project2/gilad/joycehsiao/fucci-seq/data/intensity.rds")
