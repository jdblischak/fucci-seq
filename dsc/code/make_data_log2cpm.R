#' @title Transform raw counts to log2cpm
#'
#' @param eset eset objected of final, filtered data containing counts
#'
#' @import Biobase
#'
#' @export
#'
 make_data_log2cpm < - function(eset) {
   # select endogeneous genes
   counts <- exprs(eset)[grep("ENSG", rownames(eset)), ]
   pdata <- pData(eset)
   log2cpm <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))
   return(log2cpm)
 }


