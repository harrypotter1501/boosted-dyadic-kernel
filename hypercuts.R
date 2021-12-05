# hypercuts

rm(list=ls())


library(pracma)
library(resample)
library(rdetools)


# toy data points
X = t(matrix(c(1, 3, 
             2, 5, 
             4, 1, 
             4, 7, 
             5, 2, 
             5, 4, 
             6, 1, 
             6, 6, 
             7, 4, 
             7, 5), 2, 10))
y = c(-1, -1, -1, 1, -1, 1, -1, 1, 1, 1)
plot(X)

partition = function(X, y) {
  idx = which(y == 1)
  Xp = X[idx, ]
  Xn = X[-idx, ]
  return(list(Xp=Xp, Xn=Xn))
}

normalize = function(X) {
  means = colMeans(X)
  stdevs = colStdevs(X)
  X = t(apply(X, 1, function(row) {
    row = (row - means) / stdevs
  }))
  return(list(X=X, means=means, stdevs=stdevs))
}

bdk_train = function(X, y, kernel='linear', sig=1) {
  data = normalize(X)
  Xnorm = data$X
  means = data$means
  stdevs = data$stdevs
  
  data = partition(Xnorm, y)
  Xp = data$Xp
  Xn = data$Xn
  
  if(kernel == 'rbf') {
    Kp = rbfkernel(Xnorm, sigma=sig, Y=Xp)
    Kn = rbfkernel(Xnorm, sigma=sig, Y=Xn)
  }
  else {
    Kp = Xnorm %*% t(Xp)
    Kn = Xnorm %*% t(Xn)
  }
  
  idp = which(y == 1)
  idn = which(y == -1)
  model = list(xi=X[idp[1], , drop=F], xj=X[idn[1], , drop=F], 
               kernel=kernel, sigma=sig, means=means, stdevs=stdevs, err=Inf)
  for(i in 1:ncol(Kp)) {
    ip = idp[i]
    for(j in 1:ncol(Kn)) {
      jn = idn[j]
      b = -(Kp[ip, i] - Kn[jn, j]) / 2
      pred = Kp[, i] - Kn[, j] + b
      err = norm(sigmoid(pred) - (y + 1) / 2, type='2')
      if(err < model$err) {
        xi=X[ip, , drop=F]
        xj=X[jn, , drop=F]
        model$err = err
      }
    }
  }
  
  return(model)
}

bdk_predict = function(model, X) {
  xi = model$xi
  xj = model$xj
  means = model$means
  stdevs = model$stdevs
  kernel = model$kernel
  sig = model$sigma
  
  Xnorm = t(apply(X, 1, function(row) {
    row = (row - means) / stdevs
  }))
  xi = (xi - means) / stdevs
  xj = (xj - means) / stdevs
  
  if(kernel == 'rbf') {
    Ki = rbfkernel(Xnorm, sigma=sig, Y=xi)
    Kj = rbfkernel(Xnorm, sigma=sig, Y=xj)
    b = 0
  }
  else {
    Ki = Xnorm %*% t(xi)
    Kj = Xnorm %*% t(xj)
    b = -((xi %*% t(xi) - xj %*% t(xj)) / 2)[1]
  }
  
  pred = Ki - Kj + b
  return(sigmoid(pred))
}


m = bdk_train(X, y, kernel='rbf', sig=1)

x=seq(1, 7, 0.1)
mesh = meshgrid(x)
mx = c(mesh$X)
my = c(mesh$Y)
mc = rep(0, length(mx))
Xmesh = matrix(c(mx, my), length(mx), 2)

pred = bdk_predict(m, Xmesh)

cmap=rep(3,nrow(Xmesh))
cmap[which(pred>0.5)]=5
points(Xmesh, pch=19, col=cmap)

lines(t(matrix(c(m$xi,m$xj),2,2)), lty=2)
points(t(matrix(c(m$xi,m$xj),2,2)), pch=20, col=c(4, 2))
points(X, col=y+3)

