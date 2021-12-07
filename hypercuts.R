# hypercuts

library(pracma)
library(resample)
library(rdetools)


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
      err = length(which(sign(pred) != y)) / length(y)
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

bdk_mesh = function(m, xlim, ylim, step=0.1) {
  x = seq(xlim[1], xlim[2], step)
  y = seq(ylim[1], ylim[2], step)
  mesh = meshgrid(x, y)
  mx = c(mesh$X)
  my = c(mesh$Y)
  Xmesh = matrix(c(mx, my), length(mx), 2)
  
  pred = bdk_predict(m, Xmesh)
  cmap=rep(3,nrow(Xmesh))
  cmap[which(pred>0.5)]=5
  points(Xmesh, pch=19, col=cmap)
  
  lines(t(matrix(c(m$xi,m$xj), 2, 2)), lty=2)
  points(t(matrix(c(m$xi,m$xj), 2, 2)), pch=20, col=c(4, 2))
}

ada_bdk_train = function(X, y, T, kernel='linear', sig=1) {
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
  models = list()
  for(t in 1:T) {
    #???
    for(i in 1:ncol(Kp)) {
      ip = idp[i]
      for(j in 1:ncol(Kn)) {
        jn = idn[j]
        b = -(Kp[ip, i] - Kn[jn, j]) / 2
        pred = Kp[, i] - Kn[, j] + b
        err = length(which(sign(pred) != y)) / length(y)
        if(err < model$err) {
          xi=X[ip, , drop=F]
          xj=X[jn, , drop=F]
          model$err = err
        }
      }
    }
  }
  
  return(model)
}

ada_bdk_predict = function(model, X) {
  #???
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

ada_bdk_mesh = function(m, xlim, ylim, step=0.1) {
  x = seq(xlim[1], xlim[2], step)
  y = seq(ylim[1], ylim[2], step)
  mesh = meshgrid(x, y)
  mx = c(mesh$X)
  my = c(mesh$Y)
  Xmesh = matrix(c(mx, my), length(mx), 2)
  
  pred = bdk_predict(m, Xmesh)
  cmap=rep(3,nrow(Xmesh))
  cmap[which(pred>0.5)]=5
  points(Xmesh, pch=19, col=cmap)
  
  lines(t(matrix(c(m$xi,m$xj), 2, 2)), lty=2)
  points(t(matrix(c(m$xi,m$xj), 2, 2)), pch=20, col=c(4, 2))
}

