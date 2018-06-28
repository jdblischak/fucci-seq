# evaluation



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

library(ggplot2)
ggplot(data.frame(x=10^pdata$dapi.median.log10sum.adjust),
       aes(x,y=..density..)) +
  geom_histogram(fill=NA,color="black", bins=30)+
  #  xlim(2,3.5) +
  geom_polygon(data=em.df,aes(x,y,fill=comp),color="grey50", alpha=0.5) +
  scale_fill_discrete("Component\nMeans",
                      labels=format(em.df$mu,digits=3))+
  theme_bw()
pdata$comp <- fit_mix_mclust$classification



xy_time <- lapply(1:5, function(run) {
  xy <- data.frame(
    ref_time=theta.nonvalid[part_indices[[run]]$test],
    pred_time=fits[[run]]$fit.test$cell_times_update,
    dapi=pdata$dapi.median.log10sum.adjust[match(names(theta.nonvalid[part_indices[[run]]$test]), rownames(pdata))],
    gfp=pdata$gfp.median.log10sum.adjust[match(names(theta.nonvalid[part_indices[[run]]$test]), rownames(pdata))],
    rfp=pdata$rfp.median.log10sum.adjust[match(names(theta.nonvalid[part_indices[[run]]$test]), rownames(pdata))],
    comp=pdata$comp[match(names(theta.nonvalid[part_indices[[run]]$test]), rownames(pdata))])
  return(xy)
})



for (i in 1:5) {
  xy_time[[i]]$diff_time <- pmin(
    abs(xy_time[[i]]$pred_time-xy_time[[i]]$ref_time),
    abs(xy_time[[i]]$pred_time-(2*pi-xy_time[[i]]$ref_time)))
}

plot(xy_time[[1]]$ref_time, xy_time[[1]]$dapi)




xy_time_sorted <- xy_time
# par(mfrow=c(2,3))
# for (i in 1:5) {
#   plot(xy_time_sorted[[i]]$ref_time,
#        xy_time_sorted[[i]]$dapi)
# }

cuts <- 1.86
for(run in 1:5) {
  df <- xy_time_sorted[[run]]
  for (s in 1:nrow(df)) {
    if (df$ref_time[s]>cuts) {
      df$ref_time_ord[s] <- df$ref_time[s]-cuts
    }
    if (df$ref_time[s]<=cuts) {
      df$ref_time_ord[s] <- 2*pi-(cuts-df$ref_time[s])
    }
  }
  xy_time_sorted[[run]] <- df
}

# par(mfrow=c(2,3))
# for (i in 1:5) {
#   plot(x=xy_time_sorted[[i]]$ref_time_ord,
#        y=xy_time_sorted[[i]]$pred_time)
# }
for (run in 1:5) {
  df <- xy_time_sorted[[run]]
  ii=which(df$ref_time_ord<4 & df$pred_time>2 )
  cuts <- min(df$pred_time[ii])

  for (s in 1:nrow(df)) {
    if (df$pred_time[s]>=cuts) {
      df$pred_time_ord[s] <- df$pred_time[s]-cuts
    }
    if (df$pred_time[s]<cuts) {
      df$pred_time_ord[s] <- 2*pi-(cuts-df$pred_time[s])
    }
  }
  xy_time_sorted[[run]] <- df
}


# par(mfrow=c(2,3))
# for (i in 1:5) {
#   plot(x=xy_time_sorted[[i]]$ref_time_ord,
#        y=xy_time_sorted[[i]]$pred_time_ord)
# }

for (i in 1:5) {
  xy_time_sorted[[i]]$diff_time <- pmin(
    abs(xy_time_sorted[[i]]$pred_time_ord-xy_time_sorted[[i]]$ref_time_ord),
    abs(xy_time_sorted[[i]]$pred_time_ord-(2*pi-xy_time_sorted[[i]]$ref_time_ord)))
}

par(mfrow=c(1,2))
plot(xy_time_sorted[[1]]$ref_time_ord,
     xy_time_sorted[[1]]$pred_time_ord,
     xlab="reference time", ylab="predicted time")
abline(0,1,col="blue")
plot(xy_time_sorted[[1]]$pred_time_ord,
     xy_time_sorted[[1]]$diff_time,
     xlab="predicted time", ylab="difference between times")
# plot(xy_time_sorted[[1]]$pred_time_ord,
# xy_time_sorted[[1]]$diff_time)
