# This script grabs sample metadata information from
# the imaging results


plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
            "18870_18855","18870_19101","18870_19160","19098_18511",
            "19098_18870","19098_19160","19101_18511","19101_19098",
            "19160_18870","18855_19098","19101_19160","19160_18511")

# <---- Get number of nuclei detected
dir_nuclei <- paste0("/project2/gilad/fucci-seq/intensities/")

nuclei_table <- do.call(rbind, lapply(1:length(plates), function(index) {
  plate <- plates[index]
  nnuclei <- readRDS(paste0(dir_nuclei,plate, ".nuclei.rds"))
  tmp <- data.frame(nnuclei = nnuclei,
                    well = names(nnuclei))
  tmp$plate <- plate
  return(tmp)
}))


# combine information into a data.frame
df <- data.frame(nuclei_table)

# more other data will come later...

saveRDS(df, file = "/project2/gilad/joycehsiao/fucci-seq/output/image-qc.rds")




