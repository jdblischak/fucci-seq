library(Biobase)
df <- readRDS(file="data/eset-final.rds")
pdata <- pData(df)
fdata <- fData(df)

# select endogeneous genes
counts <- exprs(df)[grep("ENSG", rownames(df)), ]

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

#macosko <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")
counts <- counts[,order(pdata$theta)]
log2cpm.all <- log2cpm.all[,order(pdata$theta)]
pdata <- pdata[order(pdata$theta),]

log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")


library(mclust)
fit_mix_mclust <- Mclust(10^pdata$dapi.median.log10sum.adjust, G=2)

# plotting
sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
x       <- seq(min(10^pdata$dapi.median.log10sum.adjust),
               max(10^pdata$dapi.median.log10sum.adjust),len=1000)
pars    <- with(fit_mix_mclust,
                data.frame(comp=factor(c(1,2)),
                           mu=parameters$mean,
                           sigma2=c(parameters$variance$sigmasq,parameters$variance$sigmasq),
                           lambda=parameters$pro))
em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
em.df$y <- with(em.df,lambda*dnorm(x,mean=mu,sd=sqrt(sigma2)))



# select external validation samples
set.seed(99)
nvalid <- round(ncol(log2cpm.quant)*.15)
ii.valid <- sample(1:ncol(log2cpm.quant), nvalid, replace = F)
ii.nonvalid <- setdiff(1:ncol(log2cpm.quant), ii.valid)

log2cpm.quant.nonvalid <- log2cpm.quant[,ii.nonvalid]
log2cpm.quant.valid <- log2cpm.quant[,ii.valid]
theta <- pdata$theta
names(theta) <- rownames(pdata)

# theta.nonvalid <- theta_moved[ii.nonvalid]
theta.nonvalid <- theta[ii.nonvalid]
theta.valid <- theta[ii.valid]

sig.genes <- readRDS("output/npreg-trendfilter-quantile.Rmd/out.stats.ordered.sig.101.rds")
expr.sig <- log2cpm.quant.nonvalid[rownames(log2cpm.quant.nonvalid) %in% rownames(sig.genes), ]


# get predicted times
# set training samples
source("peco/R/primes.R")
source("peco/R/partitionSamples.R")

parts <- partitionSamples(1:ncol(expr.sig), runs=2,
                          nsize.each = c(355, 400))
part_indices <- parts$partitions


source("peco/R/fit.trendfilter.generic.R")
source("peco/R/cycle.npreg.R")
source("code/utility.R")

Y1 <- expr.sig[,part_indices[[1]]$train]
Y2 <- expr.sig[,part_indices[[1]]$test]



# projection
theta_Y1 <- initialize_cell_times(Y1,  method.initialize.theta = "pca")
pca_Y2 <- prcomp(t(Y2), scale = TRUE)
pca_Y1_on_Y2 <- predict(pca_Y2, newdata=t(Y1))
theta_Y1_on_Y2 <- coord2rad(cbind(pca_Y1_on_Y2[,1], pca_Y1_on_Y2[,2]))
theta_Y1_on_Y2 <- as.numeric(theta_Y1_on_Y2)
plot(theta_Y1_on_Y2, theta_Y1)

# theta_rotate <- rep(0, length(theta_Y1))
# cut <- 1
# theta_rotate[which(theta_Y1_on_Y2 > cut)] <- theta_Y1_on_Y2[which(theta_Y1_on_Y2 > cut)] - cut
# theta_rotate[which(theta_Y1_on_Y2 <= cut)] <- 2*pi - theta_Y1_on_Y2[which(theta_Y1_on_Y2 <= cut)] - cut
# plot(theta_Y1_on_Y2, theta_rotate)
#
# plot(theta_Y1_on_Y2, theta_Y1)
#
# for (s in seq(0,2*pi, length.out = 100)) {
#   theta_rotate[theta_rotatetheta_Y1_on_Y2
# }

theta_Y2 <- initialize_cell_times(Y2,  method.initialize.theta = "pca")
pca_Y1 <- prcomp(t(Y1), scale = TRUE)
pca_Y2_on_Y1 <- predict(pca_Y1, newdata=t(Y2))
theta_Y2_on_Y1 <- coord2rad(cbind(pca_Y2_on_Y1[,1], pca_Y2_on_Y1[,2]))
theta_Y2_on_Y1 <- as.numeric(theta_Y2_on_Y1)
plot(theta_Y2_on_Y1, theta_Y2)



##----- method
# use it's own cell time
Y = Y1
theta = theta_Y1


initialize_grids <- function(Y, grids=100,
                            grid.method=c("pca"), ...) {

  len <- (2*pi)/(2*100)
  theta_grids <- seq(len, (2*pi)-(len), length.out=100)

  library(circular)
  if (grid.method=="pca") {
    pc_res <- prcomp(t(Y), scale = TRUE)
    grid_approx <- coord2rad(cbind(pc_res$x[,1], pc_res$x[,2]))
    grid_approx <- as.numeric(grid_approx)
    names(grid_approx) <- colnames(Y)

    theta_initial <- rep(0, length(grid_approx))
    theta_initial_ind <- rep(0, length(grid_approx))
    names(theta_initial) <- names(grid_approx)
    names(theta_initial_ind) <- names(grid_approx)

    for (i in 1:length(grid_approx)) {
    theta_initial_ind[i] <- which.min(pmin(abs(theta_grids-grid_approx[i]),
                                  abs(theta_grids-(2*pi-grid_approx[i]))))
    theta_initial <- theta_grids[theta_initial_ind]
    }
  }

  return(theta_initial)
}



method.trend = "npcirc.nw"
theta_initial <- initialize_grids(Y)
theta=theta_initial

cycle.npreg.unsupervised <- function(Y, theta,
                                 ncores=12,
                                 method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter"),
                                 maxiter=20, verbose=TRUE, tol=1,
                                 ...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial_mstep <- cycle.npreg.mstep(Y = Y_ordered,
                                     theta = theta_ordered_initial,
                                     method.trend=method.trend,
                                     ncores = ncores)

  # compute log-likelihood under initial mu and sigma
  initial_estep <- cycle.npreg.loglik(Y = initial_mstep$Y,
                                      theta = initial_mstep$theta,
                                      mu_est = initial_mstep$mu_est,
                                      sigma_est = initial_mstep$sigma_est,
                                      insample=T, funs_est=NULL)

  Y_previous <- initial_mstep$Y
  theta_previous <- initial_mstep$theta
  mu_est_previous <- initial_mstep$mu_est
  sigma_est_previous <- initial_mstep$sigma_est
  loglik_previous <- initial_estep$loglik_est
  cell_times_previous <- initial_estep$cell_times_est

  iter <- 0
  while(TRUE) {
    current_mstep <- cycle.npreg.mstep(Y = Y_previous,
                                       theta = cell_times_previous,
                                       method.trend = method.trend,
                                       ncores = ncores)
    current_loglik <- cycle.npreg.loglik(Y = current_mstep$Y,
                                theta = current_mstep$theta,
                                mu_est = current_mstep$mu_est,
                                sigma_est = current_mstep$sigma_est,
                                insample=T)

    mu_est_current <- current_mstep$mu_est
    sigma_est_current <- current_mstep$sigma_est
    funs_est_current <- current_mstep$funs
    Y_current <- current_mstep$Y
    loglik_current <- current_loglik$loglik_est
    cell_times_current <- current_loglik$cell_times_est

    if (verbose) message("log-likelihood:", loglik_current)
    eps <- loglik_current - loglik_previous

    # loop out if converged
    if (!(eps > tol & iter < maxiter)) break

    iter <- iter + 1
    loglik_previous <- loglik_current
    Y_previous <- Y_current
    mu_est_previous <- mu_est_current
    sigma_est_previous <- sigma_est_current
    cell_times_previous <- cell_times_current
    funs_est_previous <- funs_est_current
    }

  out <- list(Y_ordered=Y_current,
              loglik_est=loglik_current,
              cell_times_est=cell_times_current,
              mu_est=mu_est_current,
              sigma_est=sigma_est_current,
              funs_est=funs_est_current)
  return(out)
}

