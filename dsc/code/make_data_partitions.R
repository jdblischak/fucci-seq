
#' @param eset expression set of filtered data

#' @import Biobase
#' @export
make_data_partitions <- function(eset,
                                 training_valid=c(.85, .15),
                                 training_folds=5,
                                 type = c("individual_mixed",
                                          "individual_heldout")) {


  if (type=="individual_mixed") {

    nsamples <- ncol(eset)
    set.seed(99)
    nvalid <- round(nsamples*training_valid[2])
    ii.valid <- sample(1:nsamples, nvalid, replace = F)
    ii.train <- setdiff(1:nsamples, ii.valid)

    partitions <- partitionSamples(1:length(ii.train), runs=training_folds,
                              nsize.each = c(rep(round(length(ii.train)/5), 4),
                                             length(ii.train)-4*round(length(ii.train)/5)))

    names(partitions[[1]]) <- paste0("fold.",c(1:training_folds))
    partitions$valid <- ii.valid
    out <- list(partitions=partitions, type= "individual_mixed")
    return(out)
  }

  if (type == "individual_heldout") {

    pdata <- pData(eset)
    nsamples <- ncol(eset)
    out <- lapply(1:length(unique(pdata$chip_id)), function(j) {

      ind <- unique(pdata$chip_id)[j]
      set.seed(99)
      ii.valid <- c(1:nsamples)[which(pdata$chip_id == ind)]
      ii.train <- c(1:nsamples)[which(pdata$chip_id != ind)]

      pdata.train <- pdata[ii.train,]
      pdata.valid <- pdata[ii.valid,]

      partitions.one <- lapply(1:length(unique(pdata.train$chip_id)), function(i) {
        ind_test <- unique(pdata.train$chip_id)[i]
        test <- ii.train[which(pdata.train$chip_id==ind_test)]
        train <- ii.train[which(pdata.train$chip_id!=ind_test)]
        return(list(test=test, train=train))
      })
      names(partitions.one) <- paste0("fold.",c(1:5))
      partitions.one$valid <- ii.valid
      return(partitions.one)
    })
    names(out) <- unique(pdata$chip_id)

    out <- c(out, list(type="individual_heldout"))

    return(out)
  }

}




