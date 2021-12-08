# hypercuts

library(pracma)
library(resample)
library(rdetools)


normalize = function(X, test=FALSE, means=NULL, stdevs=NULL) {
  if(!test) {
    means = colMeans(X)
    stdevs = colStdevs(X)
  }
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
  
  if(kernel == 'rbf') {
    K = rbfkernel(Xnorm, sigma=sig)
  }
  else {
    K = Xnorm %*% t(Xnorm)
  }
  
  idp = which(y == 1)
  idn = which(y == -1)
  model = list(xi=X[idp[1], , drop=F], xj=X[idn[1], , drop=F], 
               kernel=kernel, sigma=sig, means=means, stdevs=stdevs, err=Inf)
  
  for(i in idp) {
    for(j in idn) {
      b = -(K[i, i] - K[j, j]) / 2
      pred = K[, i] - K[, j] + b
      err = length(which(sign(pred) != y)) / length(y)
      if(err < model$err) {
        model$xi = X[i, , drop=F]
        model$xj = X[j, , drop=F]
        model$err = err
      }
    }
  }
  
  return(model)
}

bdk_predict = function(model, X, likelihood=FALSE) {
  xi = model$xi
  xj = model$xj
  means = model$means
  stdevs = model$stdevs
  kernel = model$kernel
  sig = model$sigma
  
  Xnorm = normalize(X, test=TRUE, means=means, stdevs=stdevs)$X
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
  if(likelihood) {
    res = sigmoid(pred)
  }
  else {
    res = sign(pred)
  }
  return(res)
}

bdk_mesh = function(m, xlim, ylim, step=0.1) {
  x = seq(xlim[1], xlim[2], step)
  y = seq(ylim[1], ylim[2], step)
  mesh = meshgrid(x, y)
  mx = c(mesh$X)
  my = c(mesh$Y)
  Xmesh = matrix(c(mx, my), length(mx), 2)
  
  pred = bdk_predict(m, Xmesh)
  points(Xmesh, pch=19, col=pred+4)
  
  lines(t(matrix(c(m$xi,m$xj), 2, 2)), lty=2)
  points(t(matrix(c(m$xi,m$xj), 2, 2)), pch=20, col=c(4, 2))
}

w_update = function(w, a, y, pred) {
  wn = w * exp(-a * y * pred)
  wn = wn / sum(wn)
  return(wn)
}

ada_bdk_train = function(X, y, T, kernel='linear', sig=1, bs_rate=0.1) {
  data = normalize(X)
  Xnorm = data$X
  means = data$means
  stdevs = data$stdevs
  
  if(kernel == 'rbf') {
    K = rbfkernel(Xnorm, sigma=sig)
  }
  else {
    K = Xnorm %*% t(Xnorm)
  }
  
  # store results
  Xi = matrix(0, T, ncol(X))
  Xj = matrix(0, T, ncol(X))
  
  # params
  N = nrow(X)
  w = rep(1/N, N)
  a = rep(0, T)
  err = rep(0, T)
  flip = rep(1, T)
  eps = 1e-8
  
  # boost
  for(t in 1:T) {
    while(TRUE) {
      boot = sample(seq(1, N), size=floor(N * bs_rate), replace=TRUE, prob=w)
      idp = boot[which(y[boot] == 1)]
      idn = boot[which(y[boot] == -1)]
      if(length(idp) > 0 && length(idn) > 0) {
        break
      }
    }
    
    err_min = Inf
    flip_t = 1
    best_pred = rep(1, N)
    
    for(i in idp) {
      for(j in idn) {
        b = -(K[i, i] - K[j, j]) / 2
        pred = sign(K[, i] - K[, j] + b)
        err_t = sum(w[which(pred != y)])
        if(err_t > 0.5) {
          flip_t = -1
          err_t = 1 - err_t
        }
        else {
          flip_t = 1
        }
        if(err_t < err_min) {
          xi = X[i, , drop=F]
          xj = X[j, , drop=F]
          err_min = err_t
          best_pred = pred
        }
      }
    }
    
    Xi[t, ] = xi
    Xj[t, ] = xj
    
    err[t] = err_min
    flip[t] = flip_t
    a[t] = 1/2 * log((1 - err_min) / (err_min + eps))
    w = w_update(w, a[t], y, best_pred)
    
    if(err_min == 0) {
      break
    }
  }
  
  return(list(Xi=Xi[1:t, ], Xj=Xj[1:t, ], T=t, w=w, a=a[1:t], err=err[1:t], flip=flip[1:t], 
              kernel=kernel, sigma=sig, means=means, stdevs=stdevs))
}

diag_dot = function(X) {
  res = apply(X, 1, function(row) {
    return(t(row) %*% row)
  })
  return(res)
}

ada_bdk_predict = function(model, X, likelihood=FALSE) {
  means = model$means
  stdevs = model$stdevs
  kernel = model$kernel
  sig = model$sigma
  
  Xnorm = normalize(X, test=TRUE, means=means, stdevs=stdevs)$X
  Xi = normalize(model$Xi, test=TRUE, means=means, stdevs=stdevs)$X
  Xj = normalize(model$Xj, test=TRUE, means=means, stdevs=stdevs)$X
  
  if(kernel == 'rbf') {
    Ki = rbfkernel(Xnorm, sigma=sig, Y=Xi)
    Kj = rbfkernel(Xnorm, sigma=sig, Y=Xj)
    b = rep(0, nrow(X))
  }
  else {
    Ki = Xnorm %*% t(Xi)
    Kj = Xnorm %*% t(Xj)
    b = -(diag_dot(Xi) - diag_dot(Xj)) / 2
  }
  
  a = model$a
  flip = model$flip
  pred = apply(Ki - Kj + b, 1, function(row) {
    res = flip * a * row
    return(sum(res))
  })
  
  if(likelihood) {
    res = sigmoid(pred)
  }
  else {
    res = sign(pred)
  }
  
  return(res)
}

ada_bdk_mesh = function(m, xlim, ylim, step=0.1) {
  x = seq(xlim[1], xlim[2], step)
  y = seq(ylim[1], ylim[2], step)
  mesh = meshgrid(x, y)
  mx = c(mesh$X)
  my = c(mesh$Y)
  Xmesh = matrix(c(mx, my), length(mx), 2)
  
  pred = ada_bdk_predict(m, Xmesh)
  points(Xmesh, pch=19, col=pred+4)
  
  Xi = m$Xi
  Xj = m$Xj
  for(t in 1:m$T) {
    xi = Xi[t, ]
    xj = Xj[t, ]
    lines(t(matrix(c(xi, xj), 2, 2)), lty=2)
    points(t(matrix(c(xi, xj), 2, 2)), pch=20, col=c(4, 2))
  }
}

