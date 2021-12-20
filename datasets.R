###############################
# Experiments on datasets     #
###############################

rm(list=ls())

library(RCurl)

################################### WDBC #######################################
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
################################################################################
################################### WPBC #######################################
# data = read.table(
#   textConnection(
#     getURL('http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wpbc.data')
#   ),
#   sep = ',',
#   col.names=c('id_number', 'outcome', 'time', 'radius_mean',
#               'texture_mean', 'perimeter_mean', 'area_mean',
#               'smoothness_mean', 'compactness_mean',
#               'concavity_mean','concave_points_mean',
#               'symmetry_mean', 'fractal_dimension_mean',
#               'radius_se', 'texture_se', 'perimeter_se',
#               'area_se', 'smoothness_se', 'compactness_se',
#               'concavity_se', 'concave_points_se',
#               'symmetry_se', 'fractal_dimension_se',
#               'radius_worst', 'texture_worst',
#               'perimeter_worst', 'area_worst',
#               'smoothness_worst', 'compactness_worst',
#               'concavity_worst', 'concave_points_worst',
#               'symmetry_worst', 'fractal_dimension_worst',
#               'tumor_size', 'lymph_node_status')
# )
# 
# X = data.matrix(subset(
#   data,
#   select=-c(which(names(data) == 'id_number'), which(names(data) == 'outcome'))
# ))
# X <- X[complete.cases(X),]
# #X = data.matrix(data[, c('radius_mean', 'texture_mean')])
# y = rep(-1, nrow(X))
# y[which(data$outcome == 'R')] = 1 # R = recurrent, N = nonrecurrent
################################################################################
################################### WBC #######################################
# data = read.table(
#   textConnection(
#     getURL('http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data')
#   ),
#   sep = ',',
#   col.names=c('id_number', 'clump_thickness',
#               'cell_size', 'cell_shape',
#               'marginal_adhesion', 'single_epithelical_cell_size',
#               'bare_nuclei', 'bland_chromatin',
#               'normal_nucleoli',
#               'mitoses','class')
# )
# X = data.matrix(subset(
#   data,
#   select=-c(which(names(data) == 'id_number'), which(names(data) == 'class'))
# ))
# X <- X[complete.cases(X),]
# #X = data.matrix(data[, c('radius_mean', 'texture_mean')])
# y = rep(-1, nrow(X))
# y[which(data$class == '2')] = 1 #2 for benign, 4 for malignant
################################################################################
################################### WINE #######################################
# data = read.table(
#   textConnection(
#     getURL('http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data')
#   ),
#   sep = ',',
#   col.names=c('class', 'alcohol', 'malic_acid',
#               'ash', 'alcalinity',
#               'maagnesium', 'total_phenols',
#               'flavanoids', 'nonflavanoid_phenols',
#               'proanthocyanins', 'color_intensity',
#               'hue','OD280/OD315', 'proline')
# )
# X = data.matrix(subset(
#   data,
#   select=-c(which(names(data) == 'class'))
# ))
# X <- X[complete.cases(X),]
# y = rep(-1, nrow(X))
# y[which(data$class == '1')] = 1 # 1 for 1st class, 2,3 for other classes
################################################################################
################################### SONAR ######################################
# data = read.table(
#   textConnection(
#     getURL('http://archive.ics.uci.edu/ml/machine-learning-databases/undocumented/connectionist-bench/sonar/sonar.all-data')
#   ),
#   sep = ',',
# )
# X = data.matrix(subset(
#   data,
#   select=-c(61))
# )
# X <- X[complete.cases(X),]
# y = rep(-1, nrow(X))
# y[which(data[,61] == 'R')] = 1 # R for rocks, M for mines
################################################################################
################################### HEART ########################################
# data = read.table(
#   './datasets/cleveland.data', 
#   sep = ' ', 
#   fill = TRUE,
#   col.names=c('age', 'sex', 'cp', 'trestbps', 'chol', 'fbs', 'restecg', 
#               'thalach', 'exang', 'oldpeak', 'slope', 'ca', 'thal', 'num')
# )
# 
# return(data)
# 
# X = data.matrix(subset(
#   data,
#   select=-c(which(names(data) == 'num'))
# ))
# y = rep(-1, nrow(X))
# y[which(data$num == 'M')] = 1
################################################################################
# ################################## ION ########################################
# data = read.table(
#   textConnection(
#     getURL(
#       'https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data'
#     )
#   ), 
#   sep = ','
# )
# 
# X = data.matrix(subset(
#   data,
#   select=c(-2, -35)
# ))
# y = rep(-1, nrow(X))
# y[which(data[, 35] == 'g')] = 1
# 
# X[, 1:2] = as.numeric(X[, 1:2])
# ###############################################################################

sigs = c(2.5, 3, 3.5)
brs = c(0.1, 0.15, 0.2)

source('./hypercuts.R')
#params = select_params(sigs, brs, X, y)
res = ada_bdk_cv(K=5, X, y, kernel='rbf', sig=3, bs_rate=0.15, plot=TRUE)

# X2d = data.matrix(data[, c(3,4)])
# y2d = y
# plot(X2d, col=y2d+3)
# m = ada_bdk_train(X2d, y2d, X2d, y2d, T=40, sig=6)
# ada_bdk_mesh(m, c(0, 0.4), c(0, 0.5), 0.001)
# points(X2d, col=y2d+3)

