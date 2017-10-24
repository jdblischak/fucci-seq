# print number of nuclei detected
# wells that are not included in the *.nuclei.rds do
# not have any nucleus detected

plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
            "18870_18855","18870_19101","18870_19160","19098_18511",
            "19098_18870","19098_19160","19101_18511","19101_19098",
            "19160_18870")

dir_ints <- paste0("/project2/gilad/fucci-seq/intensities")

nuclei_table <- lapply(1:length(plates), function(index) {

  plate <- plates[index]

  # read in nuclei information for the plate
  fl_nuclei <- paste0(plate, ".nuclei.rds")

  # list wells with singles
  nuclei <- readRDS(file = file.path(dir_ints, fl_nuclei))
#  singles <- nuclei[which(nuclei==1)]
  data.frame(nucei_0=96-length(nuclei),
       nuclei_1=length(which(nuclei==1)),
       nuclei_2=length(which(nuclei==2)),
       nuclei_3more=length(which(nuclei>2)))
})
names(nuclei_table) <- plates
nuclei_table <- do.call(rbind, nuclei_table)

saveRDS(nuclei_table, file = "/project2/gilad/joycehsiao/fucci-seq/output/nuclei_table.rds")



# check one plate
plate <- "18855_19101"
fl_nuclei <- paste0(plate, ".nuclei.rds")

# list wells with singles
nuclei <- readRDS(file = file.path(dir_ints, fl_nuclei))

nuclei[nuclei>1]
