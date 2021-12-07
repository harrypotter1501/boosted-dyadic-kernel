
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

N = 100
X = data.matrix(data[, c('radius_mean', 'texture_mean')])[1:N, ]
y = rep(-1, nrow(X))
y[which(data$diagnosis[1:N] == 'M')] = 1

plot(X, col=y+3)

source('./hypercuts.R')
m = bdk_train(X, y, kernel='linear', sig=1)
bdk_mesh(m, c(0, 30), c(5, 40), 0.1)
points(X, col=y+3)

