source("../peco/R/utility.R")
source("../peco/R/fit.trendfilter.generic.R")

pve_wrapper <- function(results_list) {
  library(parallel)
  methods_list <- sapply(names(results_list[[1]]),
                         function(x) strsplit(x, split=".", fixed=TRUE)[[1]][2])

  res <- mclapply(1:length(results_list), function(i) {
    obj <- results_list[[i]]
    obj_res <- lapply(1:length(obj), function(j) {
      data.frame(dapi=with(obj[[j]], get.pve(dapi[order(pred_time_shift)])$pve),
                 gfp=with(obj[[j]], get.pve(gfp[order(pred_time_shift)])$pve),
                 rfp=with(obj[[j]], get.pve(rfp[order(pred_time_shift)])$pve) )
    })
    names(obj_res) <- methods_list
    out <- do.call(rbind, obj_res)
    out$fold <- names(results_list)[i]
    return(out)
  }, mc.cores=15)
  #  names(res) <- names(results_list)
  res <- do.call(rbind, res)
  res$methods <- rep(methods_list, length(results_list))
  res_mean <- aggregate(cbind(dapi, gfp, rfp)~factor(methods),
                        data=res,
                        FUN=mean)
  return(res_mean)
}
