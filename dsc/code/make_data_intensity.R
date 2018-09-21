#' @title Transform raw counts to log2cpm
#'
#' @param eset eset objected of final, filtered data containing counts
#'
#' @import Biobase
#'
#' @export
#'
 make_data_intensity < - function(eset, type=c("fucci", "fucci_dapi")) {
   # select endogeneous genes
   pdata <- pData(eset)
   if (type=="fucci") {
     mat <- data.frame(gfp=pdata$gfp.median.log10sum.adjust,
                       rfp=pdata$rfp.median.log10sum.adjust)
   }
   if (type == "fucci_dapi") {
     mat <- data.frame(gfp=pdata$gfp.median.log10sum.adjust,
                       rfp=pdata$rfp.median.log10sum.adjust,
                       dapi=pdata$dapi.median.log10sum.adjust)
   }
   return(mat)
 }


