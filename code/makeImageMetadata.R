# This script grabs sample metadata information from
# the imaging results


# <------ Image sample labels to Well labels

# read in table for matching image ID with well ID
df_label <- read.table("data/match-labels-Image2well.txt", header = TRUE)

# obtain row and column ID for the wells
# also make a table
library(stringr)
df_label$row <- str_extract(df_label$c1, "[aA-zZ]+")
df_label$column <- sprintf("%02d", as.numeric(str_extract(df_label$c1, "[0-9]+")))

df_label$image <- sprintf("%05d", df_label$image)

head(df_label)

# attach well labels to the number of nuclei detected file
plates <- c("18511_18855","18855_19101","18855_19160","18870_18511",
            "18870_18855","18870_19101","18870_19160","19098_18511",
            "19098_18870","19098_19160","19101_18511","19101_19098",
            "19160_18870","18855_19098","19101_19160","19160_18511")

# Get number of nuclei detected
dir_nuclei <- paste0("/project2/gilad/fucci-seq/intensities/")

nuclei_table <- do.call(rbind, lapply(1:length(plates), function(index) {
  plate <- plates[index]
  nnuclei <- readRDS(paste0(dir_nuclei,plate, ".nuclei.rds"))
  tmp <- data.frame(nnuclei = nnuclei,
                    image = names(nnuclei))
  tmp$plate <- plate
  return(tmp)
}))


# combine information into a data.frame
df_nuclei <- data.frame(nuclei_table)


# for some wells the image data is not included in the above nuclei data
# these are the cases that the images didn't pass the initial test for
# the presence of nuclei, I manuall add in 0 for these cases,
# need to edit the code to automatically create 0 for these cases
df_tmp <- with(df_nuclei, table(plate, nnuclei))
rowSums(df_tmp)

with(df_nuclei, table(is.na(nnuclei)))
df_nuclei <- within(df_nuclei, nnuclei[is.na(nnuclei)] <- 0)



# <------- image plate labels to experiment labels

# Read in image to experiment data label data
df_plate <- read.table("data/match-labels-Image2experiment.txt", header = TRUE)
head(df_plate)



# <------- Combine datasets
df_combo <- merge(df_nuclei, df_label, by = "image")
df_combo <- merge(df_combo, df_plate, by = "plate")

df_combo$image <- as.character(df_combo$image)
df_combo$c1 <- as.character(df_combo$c1)





# <--------- Add results of manual inspection
df_visual <- read.table("data/nuclei-visual-inspection.txt", header = TRUE,
                        stringsAsFactors = FALSE)
tmp <- sapply(df_visual$label, function(x) strsplit(x, split=".", fixed = TRUE))
df_visual$image <- sapply(tmp, "[[", 2)
df_visual$plate <- sapply(tmp, "[[", 1)

for (index in 1:nrow(df_visual)) {
  image_label <- df_visual[index,]$image
  plate_label <- df_visual[index,]$plate
  index_row <- which(df_combo$image == image_label & df_combo$plate == plate_label)
  df_combo[index_row,]$nnuclei <- 1
}


saveRDS(df_combo, file = "/project2/gilad/joycehsiao/fucci-seq/output/image-qc.rds")
write.table(df_combo, file = "/project2/gilad/joycehsiao/fucci-seq/output/image-qc.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

