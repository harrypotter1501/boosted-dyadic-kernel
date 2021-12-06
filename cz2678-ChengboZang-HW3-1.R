# hw3-1

rm(list=ls())


# prediction
ls_predict = function(X, b) {
  return(sign(X %*% b))
}

# epsilon
epsilon = function(pred, y, w) {
  return(sum(w[which(pred != y)]))
}

# base classifier
ls_classifier = function(X, y, w) {
  boot = sample(seq(1, N), replace=TRUE, prob=w)
  Xt = X[boot, ]
  yt = y[boot, ]
  b = solve(t(Xt) %*% Xt) %*% t(Xt) %*% yt
  pred = ls_predict(X, b)
  eps = epsilon(pred, y, w)
  if(eps > 0.5) {
    b = -b
    pred = -pred
    eps = epsilon(pred, y, w)
  }
  err = length(which(sign(pred) != y)) / nrow(y)
  return(list(b=b, pred=pred, eps=eps, err=err))
}

# weight update
w_update = function(w, a, y, pred) {
  wn = w * exp(-a * y * pred)
  wn = wn / sum(wn)
  return(wn)
}

# upper bound
upper_bound = function(eps) {
  res = exp(-2 * sum((1/2 - eps)^2))
  return(res)
}

# boost
adaboost = function(weights, a, X, y) {
  T = ncol(weights)
  pred = rep(0, nrow(y))
  errs = rep(0, T)
  for(t in 1:T) {
    pred = pred + a[t] * ls_predict(X, weights[, t])
    errs[t] = length(which(sign(pred) != y)) / nrow(y)
  }
  return(list(y=sign(pred), errs=errs))
}


X = data.matrix(read.csv('./hw3-data/Prob1_X.csv', header=FALSE))
y = data.matrix(read.csv('./hw3-data/Prob1_y.csv', header=FALSE))

T = 2500
N = nrow(X)

w = rep(1/N, nrow(X))
weights = matrix(0, ncol(X), T)
a = rep(0, T)
eps = rep(0, T)
bounds = rep(1, T)

w_acc = rep(0, nrow(X))

for(t in 1:T) {
  model = ls_classifier(X, y, w)
  weights[, t] = model$b
  eps[t] = model$eps
  a[t] = 1/2 * log((1 - model$eps) / model$eps)
  w = w_update(w, a[t], y, model$pred)
  bounds[t] = upper_bound(eps[1:t])
  w_acc = w_acc + w
}

boost = adaboost(weights, a, X, y)

plot(bounds, type='l', lty=2, col=3, xlab='t', ylab='error', ylim=c(0, 1))
lines(boost$errs, lty=1, col=2)
legend(1100, 0.9, legend=c('upper bound', 'training error'), lty=c(2, 1), col=c(3, 2))

barplot(c(w_acc/T), xlab='n', ylab='w')

plot(eps, type='l', col=4, xlab='t')
plot(a, type='l', col=4, xlab='t')

