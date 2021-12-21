###############################
# Dimesionality Reduction     #
###############################

rm(list=ls())

library(RCurl)

data = read.table(
  textConnection(
    getURL('https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data')
  ), 
  sep = ',', 
  col.names=c('id_number', 'diagnosis', 'radius_mean', 
              'texture_mean', 'perimeter_mean', 'area_mean', 
              'smoothness_mean', 'compactness_mean', 
              'concavity_mean','concave_points_mean', 
              'symmetry_mean', 'fractal_dimension_mean',
              'radius_se', 'texture_se', 'perimeter_se', 
              'area_se', 'smoothness_se', 'compactness_se', 
              'concavity_se', 'concave_points_se', 
              'symmetry_se', 'fractal_dimension_se', 
              'radius_worst', 'texture_worst', 
              'perimeter_worst', 'area_worst', 
              'smoothness_worst', 'compactness_worst', 
              'concavity_worst', 'concave_points_worst', 
              'symmetry_worst', 'fractal_dimension_worst')
)
PCA = function(data, n){
  library(psych)
  data_1 = principal(data, nfactors = n)
  x_2 = as.data.frame(data_1$scores)
  #cn = rep(0, n)
  cn = paste0("PC", 10:n)
  colnames(x_2) = c(cn)
  
  return (x_2)
}

SVD = function(data, n){
  x_2 = svd(data, nu = min(nrow(data), ncol(data)), nv = min(nrow(data), ncol(data)))
  for (i in 1:length(x_2$d)){
    if (x_2$d[i] < n){
      x_2$d[i] = 0
    }
    else{
      x_2$d[i] = x_2$d[i]
    }
  }
  i = which(x_2$d == 0)
  print(i)
  # for (k in 1:length(i)){
  #   x_2$u = x_2$u[,-ncol(x_2$u)]
  #   x_2$v = x_2$v[,-ncol(x_2$v)]
  # }
  # print(dim(x_2$u))
  # print(dim(x_2$v))
  #x_2$d = x_2$d[1:(i[1]-1)]
  x_2 = x_2$u%*%diag(x_2$d)%*%t(x_2$v)
  return (x_2)
}
source('./hypercuts.R')
X = data.matrix(subset(
  data,
  select=-c(which(names(data) == 'id_number'), which(names(data) == 'diagnosis'))
))
y = rep(-1, nrow(X))
y[which(data$diagnosis == 'M')] = 1
# x_2 = svd(normalize(X)$X, nu = min(nrow(data), ncol(data)), nv = min(nrow(data), ncol(data)))
# x_3 = principal(normalize(X)$X, nfactors = 8)
p = 8
s = c(0, 1, 10, 30)
a = vector(mode = "list", length = p-2+1)
res_PCA = vector(mode = "list", length = p-2+1)
res_SVD = vector(mode = "list", length = length(s))
X_index = paste0("X", 2:p)
res_index = paste0("res", 2:p)
names(a) = c(X_index)
names(res_PCA) = c(res_index)
b = vector(mode = "list", length = length(s))
names(b) = paste0("X_S_", s)
res_SVD_index = paste0("res_SVD_", s)
names(res_SVD) = c(res_SVD_index)
# for (n in 2:p){
#   a[[n-1]] = data.matrix(PCA(normalize(X)$X, n))
#   res_PCA[[n-1]] = ada_bdk_cv(K=10, a[[n-1]], y, kernel='rbf', sig=3, bs_rate=0.15, plot=FALSE)
# }
for (n in 1:length(s))
{
  b[[n]] = data.matrix(SVD(normalize(X)$X, s[n]))
  res_SVD[[n]] = ada_bdk_cv(K = 10, b[[n]], y, kernel = "rbf", sig = 3, bs_rate = 0.15, plot = FALSE)
}
#X = data.matrix(data[, c('radius_mean', 'texture_mean')])
#X = data.matrix(PCA(X, n))

sigs = c(2.5, 3, 3.5)
brs = c(0.1, 0.15, 0.2)

#params = select_params(sigs, brs, X, y)
res = ada_bdk_cv(K=10, X, y, kernel='rbf', sig=3, bs_rate=0.15, plot=TRUE)
#res_1 = ada_bdk_cv(K = 10, X_2, y, kernel = "rbf", sig = 3, bs_rate = 0.15, plot = FALSE)
#res_2 = ada_bdk_cv(K = 10, X_3, y, kernel = "rbf", sig = 3, bs_rate = 0.15, plot = FALSE)
# X2d = data.matrix(data[, c('compactness_mean', 'symmetry_mean')])
# y2d = y
# plot(X2d, col=y+3)
# m = ada_bdk_train(X2d, y2d, X2d, y2d, T=20, sig=6)
# ada_bdk_mesh(m, c(0, 0.4), c(0.05, 0.4), 0.001)
# points(X2d, col=y+3)

# plot(res$boost_errs, type='l', lty=2, col=1, xlab='rounds', ylab='err', xlim=c(1, res$T), ylim=c(0, 0.1))
# for (n in 2:p){
#   lines(res_PCA[[n-1]]$boost_errs, lty = 1, col = n)
# }
# legend(res$T-4,1,legend=c("res1", res_index), lty=c(2,rep(1, p-2+1)), col=c(1:p))
# par(c(5,5))
# plot(res$val_errs, type = "l", lty = 2, col = 1, xlab = "iterations", ylab = "error_rates", xlim=c(1, res$T), ylim=c(0, 0.15))
# for (n in 2:p){
#   lines(res_PCA[[n-1]]$val_errs, lty=1, col=n)
# }
# legend(res$T-5, 0.16, legend=c("res30",res_index), lty = c(2,rep(1, p-2+1)),col = c(1:p))
plot(res$val_errs,  type = "l", lwd = 2, lty = 2, col = 1, xlab = "iterations", ylab = "error_rates", xlim=c(1, res$T), ylim = c(0, 0.15))
for (n in 1:length(s)){
  lines(res_SVD[[n]]$val_errs, lty = 1, lwd = 2, col = n+1)
}
legend(res$T-8, 0.15, legend = c("res30", res_SVD_index), lty = c(2, 1,1,1,1), col = c(1,2,3,4,5))
# lines(res_2$upper, lty=1, col=4)
# legend((res$T-1)/2.2+1,1,legend=c('upper_raw', 'upper_PCA_8', 'upper_PCA_20'), lty=c(2, 1, 1), col=c(3, 2, 4))
