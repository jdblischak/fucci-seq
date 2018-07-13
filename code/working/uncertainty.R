results_train_top100 <- readRDS("data/results/results_train_top100.rds")

names(results_train_top100$fold.1$fit.supervised)

results_train_top100$fold.1$fit.supervised$pred_time_shift <-
  with(results_train_top100$fold.1$fit.supervised,
       rotation(ref_time, pred_time))

results_train_top100$fold.1$fit.supervised$diff_time <-
  with(results_train_top100$fold.1$fit.supervised,
       circ_dist(ref_time, pred_time_shift))

mean(results_train_top100$fold.1$fit.supervised$diff_time)/2/pi




data_training <- readRDS("data/results/data_training.rds")

sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.101.rds")

expr_sub <- data_training$log2cpm.quant.nonvalid[
  which(rownames(data_training$log2cpm.quant.nonvalid) %in% rownames(sig.genes)[1:10]),]
theta_check <- data_training$theta.nonvalid
pdata <- data_training$pdata.nonvalid

fold_indices <- readRDS("data/results/fold_indices.rds")
fold=1

Y_train <- expr_sub[,fold_indices[[fold]]$train]
theta_train <- theta_check[match(colnames(Y_train), rownames(pdata))]
names(theta_train) <- colnames(Y_train)

source("peco/R/cycle.npreg.R")
source("peco/R/fit.trendfilter.generic.R")
fit.train <- cycle.npreg.insample(Y = Y_train,
                                  theta = theta_train,
                                  ncores=15,
                                  polyorder=2,
                                  method.trend="trendfilter")
  # fitting test data
Y_test <- expr_sub[,fold_indices[[fold]]$test]
theta_test <- theta_check[match(colnames(Y_test), rownames(pdata))]
names(theta_test) <- colnames(Y_test)

fit.test <- cycle.npreg.outsample(Y_test=Y_test,
                                  sigma_est=fit.train$sigma_est,
                                  funs_est=fit.train$funs_est,
                                  method.grid = "uniform",
                                  method.trend="trendfilter",
                                  ncores=15)

fit.test$pred_time_shift <- rotation(theta_test, fit.test$cell_times_est)
fit.test$diff_time <- circ_dist(theta_test, fit.test$pred_time_shift)


credible_interval <- function(post_prob_vector, alpha=.05, verbose=F) {
  which_max <- which.max(post_prob_vector)
  if (which_max==length(post_prob_vector)) {
    forward_grid_prob <- post_prob_vector[1:(which_max-1)]
    reverse_grid_prob <- post_prob_vector[rev(1:(which_max-1))]
  } else if (which_max==1) {
    forward_grid_prob <- post_prob_vector[2:length(post_prob_vector)]
    reverse_grid_prob <- post_prob_vector[rev(2:length(post_prob_vector))]
  } else {
    forward_grid_prob <- post_prob_vector[c(c((which_max+1):length(post_prob_vector)),
                                            c(1:(which_max-1)))]
    reverse_grid_prob <- post_prob_vector[c(rev(1:(which_max-1)),
                                            rev((which_max+1):length(post_prob_vector)))]
  }

  sum_post_prob_initial <- post_prob_vector[which_max]
  i=1

  while(i<length(forward_grid)) {
    #print(i)
    sum_post_prob_initial <- sum_post_prob_initial + forward_grid_prob[i] + reverse_grid_prob[i]
    sum_post_prob_update <- sum_post_prob_initial

    if (verbose) message("cum. prob:", sum_post_prob_update)
    grid <- i

    if (sum_post_prob_update > (1-alpha)) break
    i=i+1
    sum_post_prob_initial <- sum_post_prob_update

    }

    return(list(iter=i,
                sum_post_prob=sum_post_prob_update,
                forward_boundary=as.numeric(names(forward_grid_prob[grid])),
                reverse_boundary=as.numeric(names(reverse_grid_prob[grid]))))
}



grids <- as.numeric(colnames(fit.test$prob_per_cell_by_celltimes))
ncells <- ncol(fit.test$Y)
dd <- rep(NULL, ncells)
for (j in 1:ncells) {
  #print(j)
  #cell_ind <- order(fit.test$diff_time, decreasing = T)[j]
  post_prob_vector <- fit.test$prob_per_cell_by_celltimes[j,]
  rr=credible_interval(post_prob_vector, verbose=F)
  which_max <- which.max(post_prob_vector)
  forward_grids <- c(c((which_max+1):length(post_prob_vector)),
                     c(1:(which_max-1)))
  reverse_grids <- c(rev(1:(which_max-1)),
                     rev((which_max+1):length(post_prob_vector)))
  cols <- rep("gray80", 100)
  cols[c(reverse_grids[1:rr$iter])] <- "lightblue"
  cols[c(forward_grids[1:rr$iter])] <- "lightblue"
  cols[which_max] <- "lightblue"
  grids_outside <- grids[which(cols=="gray80")]
  # TRUE inside credible interval
  # FALSE otuside credible interval
  dd[j] <- !(theta_test[j] >= min(grids_outside) &
      theta_test[j] <= max(grids_outside))
}


sum(dd)
which(dd==F)
par(mfrow=c(1,1))
plot(theta_test, fit.test$pred_time_shift,
     col=c("red", "gray70")[dd+1], pch=16,
     main = "Fucci time outside/inside (red/gray) credible itnerval")
mean(dd)

par(mfrow=c(2,5))
for (j in 1:10) {
  cell_ind <- which(dd==F)[j]
  post_prob_vector <- fit.test$prob_per_cell_by_celltimes[cell_ind,]
  rr=credible_interval(post_prob_vector, verbose=F)
  which_max <- which.max(post_prob_vector)
  forward_grids <- c(c((which_max+1):length(post_prob_vector)),
                     c(1:(which_max-1)))
  reverse_grids <- c(rev(1:(which_max-1)),
                     rev((which_max+1):length(post_prob_vector)))
  cols <- rep("gray80", 100)
  cols[c(reverse_grids[1:rr$iter])] <- "lightblue"
  cols[c(forward_grids[1:rr$iter])] <- "lightblue"
  cols[which_max] <- "lightblue"
  # grids_outside <- grids[which(cols=="gray80")]
  # if (theta_test[cell_ind] >= min(grids_outside) &
  #     theta_test[cell_ind] <= max(grids_outside))
  plot(fit.test$prob_per_cell_by_celltimes[cell_ind,], type="h",
       xaxt="n", ylab="Posterior prob", col=cols, axes=F,
       xlab="Grid points", lwd=3)
  axis(1, at=c(1,seq(20, 100, by=20)),
       labels=round(grids[c(1,seq(20, 100, by=20))],3))
  axis(2)
  abline(v=which.max(fit.test$prob_per_cell_by_celltimes[cell_ind,]),
         col="darkblue", lty=1, lwd=3)
  abline(v=which.min(circ_dist(grids, theta_test[cell_ind])),
         col="red", lty=1, lwd=3)
  title(paste(round(fit.test$diff_time[cell_ind]/2/pi,2)))
}



for (j in 1:5) {
  cell_ind <- order(fit.test$diff_time, decreasing = F)[j]
  post_prob_vector <- fit.test$prob_per_cell_by_celltimes[cell_ind,]
  rr=credible_interval(post_prob_vector, verbose=F)
  which_max <- which.max(post_prob_vector)
  forward_grids <- c(c((which_max+1):length(post_prob_vector)),
                     c(1:(which_max-1)))
  reverse_grids <- c(rev(1:(which_max-1)),
                     rev((which_max+1):length(post_prob_vector)))
  cols <- rep("gray80", 100)
  cols[c(reverse_grids[1:rr$iter])] <- "lightblue"
  cols[c(forward_grids[1:rr$iter])] <- "lightblue"
  plot(fit.test$prob_per_cell_by_celltimes[cell_ind,], type="h",
       xaxt="n", ylab="Posterior prob", col=cols, axes=F,
       xlab="Grid points", lwd=3)
  axis(1, at=c(1,seq(20, 100, by=20)),
       labels=round(grids[c(1,seq(20, 100, by=20))],3))
  axis(2)
  abline(v=which.max(fit.test$prob_per_cell_by_celltimes[cell_ind,]),
         col="darkblue", lty=1, lwd=3)
  abline(v=which.min(circ_dist(grids, theta_test[cell_ind])),
         col="red", lty=1, lwd=3)
  title(paste(round(fit.test$diff_time[cell_ind]/2/pi,2)))
}


