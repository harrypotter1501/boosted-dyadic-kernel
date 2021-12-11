
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

M = 1
N = 100
X = data.matrix(subset(
  data,
  select=-c(which(names(data) == 'id_number'), which(names(data) == 'diagnosis'))
))
#X = data.matrix(data[, c('radius_mean', 'texture_mean')])
y = rep(-1, nrow(X))
y[which(data$diagnosis == 'M')] = 1

#plot(X, col=y+3)

source('./hypercuts.R')
m = ada_bdk_train(X, y, T=100, kernel='rbf', sig=2, bs_rate=0.1)
err = length(which(ada_bdk_predict(m, X) != y)) / length(y)
# ada_bdk_mesh(m, c(0, 30), c(5, 40), 0.1)
# points(X, col=y+3)

plot(m$upper, type='l', lty=2, col=3, xlab='boost rounds', ylab='err', xlim=c(1, m$T), ylim=c(0, 1))
lines(m$boost_errs, lty=1, col=2)
legend((m$T-1)/2.2+1, 1.0, legend=c('upper bound', 'training error'), lty=c(2, 1), col=c(3, 2))

