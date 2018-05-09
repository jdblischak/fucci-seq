library(genlasso)
# library(NPCirc)
# library(smashr)
# yy <- yy.train
# pos.yy <- pos.yy.train

fit.trendfilter.generic <- function(yy, pos.yy=c(1:length(yy)), polyorder=2) {

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  pos.rep <- rep(pos.yy,3)
  include <- rep(c(FALSE, TRUE, FALSE), each = n)

  # trendfilter
  fit.trend <- trendfilter(yy.rep,
                           ord=polyorder, approx = F, maxsteps = 1000)
  cv.trend <- cv.trendfilter(fit.trend)
  which.lambda <- cv.trend$i.1se
  yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                           df=fit.trend$df[which.lambda])$fit
  trend.yy <- yy.trend.pred[include]
  pve <- 1-var(yy-trend.yy)/var(yy)

  # trend.mad.pred <- mad(yy.nonzero.rep[include]-yy.trend.pred[include], constant = 1)
  # trend.mad.constant <- mad(yy.nonzero.rep[include]-mean(yy.nonzero.rep[include]), constant = 1)
  # trend.mad.ratio <- trend.mad.pred/trend.mad.constant

  return(list(trend.yy=trend.yy,
              trend.pos=pos.rep[include],
              #trend.mad.ratio=trend.mad.ratio,
              trend.pve=pve))
}

fit.trendfilter <- function(yy, pos.yy=c(1:length(yy)), polyorder=2) {
  nonzero <- which(yy!=0)
  yy.nonzero <- yy[nonzero]
  pos.nonzero <- pos.yy[nonzero]

  n <- length(yy.nonzero)
  yy.nonzero.rep <- rep(yy.nonzero,3)
#  theta.nonzero.rep <- rep(theta.nonzero,3)
  pos.nonzero.rep <- rep(pos.nonzero,3)
  include <- rep(c(FALSE, TRUE, FALSE), each = n)

  # trendfilter
#  fit.trend <- trendfilter(yy.nonzero, ord=2, pos = pos.nonzero, maxsteps = 1000)
  fit.trend <- trendfilter(yy.nonzero.rep,
                           ord=polyorder, approx = F, maxsteps = 1000)
  cv.trend <- cv.trendfilter(fit.trend)
  which.lambda <- cv.trend$i.1se
  yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                           df=fit.trend$df[which.lambda])$fit
  trend.yy <- yy.trend.pred[include]
  pve <- 1-var(yy-trend.yy)/var(yy)

  # trend.mad.pred <- mad(yy.nonzero.rep[include]-yy.trend.pred[include], constant = 1)
  # trend.mad.constant <- mad(yy.nonzero.rep[include]-mean(yy.nonzero.rep[include]), constant = 1)
  # trend.mad.ratio <- trend.mad.pred/trend.mad.constant

  return(list(trend.yy=trend.yy,
              trend.pos=pos.nonzero.rep[include],
              #trend.mad.ratio=trend.mad.ratio,
              trend.pve=pve))
}




fit.trendfilter.includezero <- function(yy, pos.yy=c(1:length(yy)), polyorder=2) {
#  nonzero <- which(yy!=0)
  nonzero <- c(1:length(yy))
  yy.nonzero <- yy[nonzero]
  pos.nonzero <- pos.yy[nonzero]

  n <- length(yy.nonzero)
  yy.nonzero.rep <- rep(yy.nonzero,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  pos.nonzero.rep <- rep(pos.nonzero,3)
  include <- rep(c(FALSE, TRUE, FALSE), each = n)

  # trendfilter
  #  fit.trend <- trendfilter(yy.nonzero, ord=2, pos = pos.nonzero, maxsteps = 1000)
  fit.trend <- trendfilter(yy.nonzero.rep,
                           ord=polyorder, approx = F, maxsteps = 2000)
  cv.trend <- cv.trendfilter(fit.trend)
  which.lambda <- cv.trend$i.1se
  yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                           df=fit.trend$df[which.lambda])$fit

  trend.mad.pred <- mad(yy.nonzero.rep[include]-yy.trend.pred[include], constant = 1)
  trend.mad.constant <- mad(yy.nonzero.rep[include]-mean(yy.nonzero.rep[include]), constant = 1)
  trend.mad.ratio <- trend.mad.pred/trend.mad.constant

  return(list(trend.yy=yy.trend.pred[include],
              trend.pos=pos.nonzero.rep[include],
              trend.mad.ratio=trend.mad.ratio,
              trend.cv=cv.trend$i.1se,
              trend.df=fit.trend$df[which.lambda],
              trend.lambda=cv.trend$lambda.1se))
}


fit.kernel <- function(yy, theta) {
  yy.nonzero <- yy[yy!=0]
  theta.nonzero <- theta[yy!=0]

  # local linear kernel estimate
  fit.npcirc.nw <- kern.reg.circ.lin(x = theta.nonzero,
                                     y = yy.nonzero, method="NW")
  npnw.lfun <- approxfun(y=fit.npcirc.nw$y,x=fit.npcirc.nw$x, rule=2)
  npnw.pred.y <- npnw.lfun(theta.nonzero)
  npnw.mad.pred <- mad(yy.nonzero-npnw.pred.y, constant = 1, na.rm=TRUE)
  npnw.mad.constant <- mad(yy.nonzero-mean(yy.nonzero), constant = 1)
  npnw.mad.ratio <- npnw.mad.pred/npnw.mad.constant
  #  npll.pve <- 1-var(yy.train.nonzero-npll.pred.y.train)/var(yy.train.nonzero)

  return(list(npnw.yy=npnw.pred.y,
              npnw.pos=c(1:length(yy))[yy!=0],
              npnw.theta=theta.nonzero,
              npnw.mad.ratio=npnw.mad.ratio))
}



fit.smash <- function(yy, theta) {
  yy <- yy[order(theta)]
  is.zeros <- which(yy == 0)
  yy.impute <- yy
  yy.impute[yy==0] <- mean(yy[yy!=0])

  fit.smash <- smash(x = yy.impute, model="gaus")

  # plot(c(1:length(yy)), yy)
  # points(c(1:length(fit.smash)), fit.smash, col = "royalblue", pch=16, cex=.7)
  smash.mad.pred <- mad(yy.impute-fit.smash, constant = 1)
  smash.mad.constant <- mad(yy.impute-mean(yy.impute), constant = 1)
  smash.mad.ratio <- smash.mad.pred/smash.mad.constant

  return(list(smash.yy=fit.smash,
              smash.pos=c(1:length(yy)),
              smash.mad.ratio=smash.mad.ratio))
}

