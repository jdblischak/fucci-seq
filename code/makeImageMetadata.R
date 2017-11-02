# This script grabs sample metadata information from
# the imaging results


plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
            "18870_18855","18870_19101","18870_19160","19098_18511",
            "19098_18870","19098_19160","19101_18511","19101_19098",
            "19160_18870","18855_19098","19101_19160","19160_18511")

# <---- Add a column of the C1 well labels
df_label <- data.frame(well=c(1:96))
df_label$well_c1 <- df_label$well
df_label <- within(df_label, {
  well_c1[well == 3] <- "A1"
  well_c1[well == 2] <- "A2"
  well_c1[well == 1] <- "A3"
  well_c1[well == 49] <- "A4"
  well_c1[well == 50] <- "A5"
  well_c1[well == 6] <- "A6"
  well_c1[well == 5] <- "A7"
  well_c1[well == 4] <- "A8"
  well_c1[well == 52] <- "A9"
  well_c1[well == 53] <- "A10"
  well_c1[well == 54] <- "A11"
  well_c1[well == 12] <- "A12"
  well_c1[well == 9] <- "B1"
  well_c1[well == 8] <- "B2"
  well_c1[well == 7] <- "B3"
  well_c1[well == 55] <- "B4"
  well_c1[well == 56] <- "B5"
  well_c1[well == 57] <- "B6"
  well_c1[well == 12] <- "B7"
  well_c1[well == 11] <- "B8"
  well_c1[well == 10] <- "B9"
  well_c1[well == 58] <- "B10"
  well_c1[well == 59] <- "B11"
  well_c1[well == 60] <- "B12"
  well_c1[well == 15] <- "C1"
  well_c1[well == 14] <- "C2"
  well_c1[well == 13] <- "C3"
  well_c1[well == 61] <- "C4"
  well_c1[well == 62] <- "C5"
  well_c1[well == 63] <- "C6"
  well_c1[well == 18] <- "C7"
  well_c1[well == 17] <- "C8"
  well_c1[well == 16] <- "C9"
  well_c1[well == 64] <- "C10"
  well_c1[well == 65] <- "C11"
  well_c1[well == 66] <- "C12"
  well_c1[well == 21] <- "D1"
  well_c1[well == 20] <- "D2"
  well_c1[well == 19] <- "D3"
  well_c1[well == 67] <- "D4"
  well_c1[well == 68] <- "D5"
  well_c1[well == 69] <- "D6"
  well_c1[well == 24] <- "D7"
  well_c1[well == 23] <- "D8"
  well_c1[well == 22] <- "D9"
  well_c1[well == 70] <- "D10"
  well_c1[well == 71] <- "D11"
  well_c1[well == 72] <- "D12"
  well_c1[well == 25] <- "E1"
  well_c1[well == 26] <- "E2"
  well_c1[well == 27] <- "E3"
  well_c1[well == 75] <- "E4"
  well_c1[well == 74] <- "E5"
  well_c1[well == 73] <- "E6"
  well_c1[well == 28] <- "E7"
  well_c1[well == 29] <- "E8"
  well_c1[well == 30] <- "E9"
  well_c1[well == 78] <- "E10"
  well_c1[well == 77] <- "E11"
  well_c1[well == 76] <- "E12"
  well_c1[well == 31] <- "F1"
  well_c1[well == 32] <- "F2"
  well_c1[well == 33] <- "F3"
  well_c1[well == 81] <- "F4"
  well_c1[well == 80] <- "F5"
  well_c1[well == 79] <- "F6"
  well_c1[well == 34] <- "F7"
  well_c1[well == 35] <- "F8"
  well_c1[well == 36] <- "F9"
  well_c1[well == 84] <- "F10"
  well_c1[well == 83] <- "F11"
  well_c1[well == 82] <- "F12"
  well_c1[well == 37] <- "G1"
  well_c1[well == 38] <- "G2"
  well_c1[well == 39] <- "G3"
  well_c1[well == 87] <- "G4"
  well_c1[well == 86] <- "G5"
  well_c1[well == 85] <- "G6"
  well_c1[well == 40] <- "G7"
  well_c1[well == 41] <- "G8"
  well_c1[well == 42] <- "G9"
  well_c1[well == 90] <- "G10"
  well_c1[well == 89] <- "G11"
  well_c1[well == 88] <- "G12"
  well_c1[well == 43] <- "H1"
  well_c1[well == 44] <- "H2"
  well_c1[well == 45] <- "H3"
  well_c1[well == 93] <- "H4"
  well_c1[well == 92] <- "H5"
  well_c1[well == 91] <- "H6"
  well_c1[well == 46] <- "H7"
  well_c1[well == 47] <- "H8"
  well_c1[well == 48] <- "H9"
  well_c1[well == 96] <- "H10"
  well_c1[well == 95] <- "H11"
  well_c1[well == 94] <- "H12"
})

df_label$well <- sprintf("%05d", df_label$well)

library(stringr)
df_label$well_c1_row <- str_extract(df_label$well_c1, "[aA-zZ]+")
df_label$well_c1_row <- str_extract(df_label$well_c1, "[aA-zZ]+")
df_label$well_c1_column <- sprintf("%02d", as.numeric(str_extract(df_label$well_c1, "[0-9]+")))

df_label$well_c1 <- paste0(df_label$well_c1_row, df_label$well_c1_column)



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
df_nuclei <- data.frame(nuclei_table)



#<---  combine datasets
df_combo <- merge(df_nuclei, df_label, by = "well")
saveRDS(df_combo, file = "/project2/gilad/joycehsiao/fucci-seq/output/image-qc.rds")



