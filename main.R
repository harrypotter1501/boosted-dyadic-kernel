
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

X = data.matrix(subset(
  data,
  select=-c(which(names(data) == 'id_number'), which(names(data) == 'diagnosis'))
))
#X = data.matrix(data[, c('radius_mean', 'texture_mean')])
y = rep(-1, nrow(X))
y[which(data$diagnosis == 'M')] = 1


sigs = c(2.5, 3, 3.5)
brs = c(0.1, 0.15, 0.2)

source('./hypercuts.R')
#params = select_params(sigs, brs, X, y)
res = ada_bdk_cv(K=5, X, y, kernel='rbf', sig=3, bs_rate=0.15, plot=TRUE)

# X2d = data.matrix(data[, c('compactness_mean', 'symmetry_mean')])
# y2d = y
# plot(X2d, col=y+3)
# m = ada_bdk_train(X2d, y2d, X2d, y2d, T=20, sig=6)
# ada_bdk_mesh(m, c(0, 0.4), c(0.05, 0.4), 0.001)
# points(X2d, col=y+3)

