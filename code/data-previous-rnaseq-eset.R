# Leng et al. 2015
# Download data from my GitHub repo
install_github("jhsiao999/singleCellRNASeqHumanLengESC")
library("singleCellRNASeqHumanLengESC")
data("HumanLengESC")
saveRDS(HumanLengESC, file = "data/rnaseq-previous-studies/HumanLengESC.rds")

