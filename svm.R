###############################
# SVM with RBF Kernel         #
###############################

rm(list=ls())

library(RCurl)
library(e1071)
library(caTools)

data = read.table(
  textConnection(
    getURL(
      'https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data'
    )
  ), 
  sep = ','
)
X = data.matrix(subset(
  data,
  select=-c(2,35))
)
X <- X[complete.cases(X),]
# Select two attributes
#X = data.matrix(data[, c('radius_mean', 'texture_mean')])
y = rep(-1, nrow(X))
y[which(data[,35] == 'g')] = 1 # g for good

plot(X, col=y+3)

# Data Set Organization
split = sample.split(seq(1,nrow(X)), SplitRatio = 0.8)
X_train = subset(X, split == TRUE)
X_val = subset(X, split == FALSE)
y_train = subset(y, split == TRUE)
y_val = subset(y, split == FALSE)

svmfit = svm(y_train ~ ., data = X_train, type = 'C-classification', kernel = "radial", cost = 10, scale = FALSE)
print(svmfit)

y_pred = predict(svmfit, X_val)
confusion = table(y_val,y_pred)
confusion
err_test = length(which(y_pred != y_val)) / length(y_val)
print(sprintf('Test Accuracy = %f', (1 - err_test)*100))

# Plot
#plot(X_val, col=y_val+3)

# 10-Fold Cross-Validation
#svm_cv = tune.svm(X, matrix(y, nrow(X), 1), type = 'C-classification', kernel = "radial", cost = c(0.001,0.01,0.1,10,100), scale = FALSE)
#X_train.cv <- X_train[complete.cases(X_train),]
#X_train.cv <- as.data.frame(X_train.cv)
#y_train.cv <- y_train[complete.cases(y_train)]
#svm.cv = tune(svm, y_train.cv ~ ., data = X_train.cv, type = 'C-classification', kernel = "radial", cost = c(0.001,0.01,0.1,1,5,10,100,1000), scale = FALSE)
