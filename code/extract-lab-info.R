library("dplyr")
x <- read.delim("../data/fucci_qc.txt", stringsAsFactors = FALSE,
                colClasses = c("character", "character", "numeric",
                               "numeric", "character", "character",
                               "character", "character", "character"))

stopifnot(x$ERCC == "1:50000")

out <- x %>%
  mutate(sample = paste(experiment, well, sep = "-"),
         ERCC = "50x dilution") %>%
  select(sample, experiment:image_label) %>%
  arrange(sample)

dir.create("../data/lab-info", showWarnings = FALSE)

for (e in unique(out$experiment)) {
  print(e)
  fname <- paste0("../data/lab-info/", e, ".txt")
  d <- out %>% filter(experiment == e)
  write.table(d, file = fname, quote = FALSE, sep = "\t", row.names = FALSE)
}
