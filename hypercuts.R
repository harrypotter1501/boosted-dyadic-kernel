###############################
# Model Definition and Utils  #
###############################

library(pracma)
library(resample)
library(rdetools)
library(caTools)


# Normalize each dimension
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

# Adaboost training error upper bound
upper_bound = function(eps) {
  res = exp(-2 * sum((1/2 - eps)^2))
  return(res)
}

# Adaboost update weights
w_update = function(w, a, y, pred) {
  wn = w * exp(-a * y * pred)
  wn = wn / sum(wn)
  return(wn)
}

# Train
ada_bdk_train = function(X, y, X_val, y_val, T=30, kernel='rbf', sig=3, bs_rate=0.15) {
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
  eps = 1e-16
  
  # model for every round
  model = list()
  
  # errors
  boost_errs = rep(0, T)
  val_errs = rep(0, T)
  upper = rep(0, T)
  
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
    best_pred = rep(1, N)
    
    for(i in idp) {
      for(j in idn) {
        b = -(K[i, i] - K[j, j]) / 2
        pred = K[, i] - K[, j] + b
        err_t = sum(w[which(sign(pred) != y)])
        
        if(err_t < err_min) {
          xi = X[i, , drop=F]
          xj = X[j, , drop=F]
          best_pred = pred
          flip[t] = 1
          err_min = err_t
        }
        else if(1 - err_t < err_min) {
          xi = X[i, , drop=F]
          xj = X[j, , drop=F]
          err_min = 1 - err_t
          best_pred = -pred
          flip[t] = -1
        }
      }
    }
    
    Xi[t, ] = xi
    Xj[t, ] = xj
    
    err[t] = err_min
    a[t] = 1/2 * log((1 - err_min) / (err_min + eps))
    w = w_update(w, a[t], y, sign(best_pred))
    
    # test
    model = list(Xi=Xi[1:t, ,drop=F], Xj=Xj[1:t, , drop=F], T=t, w=w, a=a[1:t], err=err[1:t], flip=flip[1:t], 
                 kernel=kernel, sigma=sig, means=means, stdevs=stdevs)
    
    be = length(which(ada_bdk_predict(model, X) != y)) / length(y)
    boost_errs[t] = be
    upper[t] = upper_bound(err[1:t])
    
    # val
    val_be = length(which(ada_bdk_predict(model, X_val) != y_val)) / length(y_val)
    val_errs[t] = val_be
    
    model = append(model, list(boost_errs=boost_errs[1:t], val_errs=val_errs[1:t], upper=upper[1:t]))
    #print(sprintf('Boosting %d rounds... Current w_err = %f, train_acc = %f, val_acc = %f', t, err_min, 1 - be, 1 - val_be))
    
    if(be == 0) {
      break
    }
  }
  
  return(model)
}

# Calcuate dot products along the diagonals
diag_dot = function(X) {
  res = apply(X, 1, function(row) {
    return(t(row) %*% row)
  })
  return(res)
}

# Prediction
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
    b = 0
  }
  else {
    Ki = Xnorm %*% t(Xi)
    Kj = Xnorm %*% t(Xj)
    b = -(diag_dot(Ki) - diag_dot(Kj)) / 2
  }
  
  a = model$a
  flip = model$flip
  
  pred = apply(Ki - Kj + b, 1, function(row) {
    return(sum(flip * a * sign(row)))
  })
  
  if(likelihood) {
    res = sigmoid(pred)
  }
  else {
    res = sign(pred)
  }
  
  return(res)
}

# Plot 2D example
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

# Crossvalidation
ada_bdk_cv = function(K=10, X, y, T=100, kernel='rbf', sig=3, bs_rate=0.15, plot=TRUE) {
  N = nrow(X)
  ids = sample(seq(1, N), N)
  val_size = floor(N / K)
  
  upper = rep(0, T)
  boost_errs = rep(0, T)
  val_errs = rep(0, T)
  t = T
  
  for(k in 1:K) {
    val_id = seq((k - 1) * val_size, k * val_size)
    X_val = X[val_id, ]
    y_val = y[val_id]
    X_train = X[-val_id, ]
    y_train = y[-val_id]
    
    m = ada_bdk_train(X_train, y_train, X_val, y_val, T=T, kernel=kernel, sig=sig, bs_rate=bs_rate)
    
    if(m$T < t) {
      t = m$T
    }
    upper[1:t] = upper[1:t] + m$upper[1:t] / K
    boost_errs[1:t] = boost_errs[1:t] + m$boost_errs[1:t] / K
    val_errs[1:t] = val_errs[1:t] + m$val_errs[1:t] / K
    
    if(plot) {
      print(sprintf('Validating %d folds...', k))
    }
  }
  
  if(plot) {
    plot(upper[1:t], type='l', lty=2, col=3, xlab='rounds', ylab='err', xlim=c(1, t), ylim=c(0, 1))
    lines(boost_errs[1:t], lty=1, col=2)
    lines(val_errs[1:t], lty=1, col=4)
    legend((t-1)/2.2+1, 1.0, legend=c('upper_bound', 'train_error', 'val_error'), lty=c(2, 1, 1), col=c(3, 2, 4))
  }
  
  return(list(T=t, boost_errs=boost_errs[1:t], val_errs=val_errs[1:t]))
}


# CV select parameter
select_params = function(sigs, bs_rates, X, y, K=5, T=30, kernel='rbf') {
  val_err_min = Inf
  sig_cv = sigs[1]
  bs_cv = bs_rates[1]
  res_cv = NULL
  t_cv = 0
  
  plot(-1, -1, xlab='rounds', ylab='err', xlim=c(1, 30), ylim=c(0, 0.15))
  
  for(s in sigs) {
    for(r in bs_rates) {
      print(sprintf('Validating sigma = %f & bs_rate = %f...', s, r))
      
      t = Sys.time()
      res = ada_bdk_cv(K, X, y, T=T, kernel=kernel, sig=s, bs_rate=r, plot=FALSE)
      te = difftime(Sys.time(), t, units='secs')[[1]]
      val_err = res$val_errs[res$T]
      
      if(val_err < val_err_min) {
        sig_cv = s
        bs_cv = r
        res_cv = res
        t_cv = te
        val_err_min = val_err
      }
      
      print(sprintf('Finished within %d rounds with acc = %f, time elapsed %f', res$T, 1 - val_err, te))
      
      lines(res$boost_errs, lty=1, col=2)
      lines(res$val_errs, lty=1, col=4)
    }
  }
  
  legend((30-1)/1.8+1, 0.15, legend=c('train_err', 'val_err'), lty=c(1, 1), col=c(2, 4))
  
  return(list(T=res_cv$T, sig=sig_cv, bs_rate=bs_cv, t=t_cv))
}

