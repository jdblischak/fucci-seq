# Description:
#
#   This script imports the nuclei information rds from the image analysis.
# Modify the results based on visual inspection. Save the final data to
# image-qc-final.rds.

image_qc <- readRDS("output/image-qc.rds")
image_qc$plate <- as.character(image_qc$plate)
image_qc$image <- as.character(image_qc$image)


library(data.table)
multiples <- fread("data/nuclei-visual-inspection-multiple.txt",
                   header = FALSE, stringsAsFactors = FALSE)
zeros <- fread("data/nuclei-visual-inspection-zero.txt",
               header = FALSE, stringsAsFactors = FALSE)

multiples_list <- lapply(1:nrow(multiples), function(index) {
  strsplit(multiples$V1[index], split=".", fixed=TRUE)[[1]]
})

zeros_list <- lapply(1:nrow(zeros), function(index) {
  strsplit(zeros$V1[index], split=".", fixed=TRUE)[[1]]
})

for (index in 1:length(multiples_list)) {
  image_qc[with(image_qc,
      which(plate==multiples_list[[1]][1] & image == multiples_list[[1]][2])), ] <- 1
}

for (index in 1:length(zeros_list)) {
  image_qc[with(image_qc,
                which(plate==zeros_list[[1]][1] & image == zeros_list[[1]][2])), ] <- 1
}

saveRDS(image_qc, file = "data/image-qc-final.rds")
