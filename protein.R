cat("PROTEIN")
set.seed (1000)
library(glmnet)
library(plsgenomics)
library(pls)
library(e1071)

load("path/protein.RData")
Y_FULL = data.train[,"Y"]
X_FULL = data.train[,0:2000]

#CORRELATION ON HEATMAP
correlation_full <- cor(na.omit(X_FULL))
png('protein_heatmap.png')
matrix.heatmap(correlation_full)
dev.off()


FULL = data.train[sample(nrow(data.train)),]
folds = cut(seq(1,nrow(FULL)),breaks=10,labels=FALSE)

ridge_R2_mean=0
ridge_MSE_mean=0
lasso_R2_mean=0
lasso_MSE_mean=0


#Self made cross validation loop
for(i in 1:10){
  cat("\nValidation number: ", i)

  testIndexes <- which(folds==i,arr.ind=TRUE)

  test <- FULL[testIndexes, ]
  X_test = test[,0:2000]
  Y_test = test[,"Y"]

  train <- FULL[-testIndexes, ]
  Y = train[,"Y"]
  X = train[,0:2000]

  #LASSO REGRESSION METHOD
  lasso_lambda = cv.glmnet(y=Y,x=X,alpha = 1)
  lasso_model = glmnet(X,Y,alpha=1,lambda=lasso_lambda$lambda.min)
  lasso_predict = predict(lasso_model,s=lasso_lambda$lambda.min,newx=X_test)
  lasso_coef = which(as.data.frame(as.matrix(coef(lasso_model))) > 0, arr.ind = TRUE)
  lasso_MSE=sum((Y_test-lasso_predict)^2)/length(Y_test)
  lasso_R2= 1 - (sum((Y_test-lasso_predict )^2)/sum((Y_test-mean(Y_test))^2))
  cat('\nLasso R2: ',lasso_R2)
  cat("\nLasso MSE: ",lasso_MSE)
  
  #RIDGE REGRESSION METHOD
  ridge_lambda = cv.glmnet(X, Y, alpha = 0)
  ridge_model = glmnet(X, Y, alpha=0, lambda = ridge_lambda$lambda.min)
  ridge_predict = predict(ridge_model,s=ridge_lambda$lambda.min,newx=X_test)
  ridge_coef = as.matrix(coef(ridge_model))
  ridge_MSE=sum((Y_test-ridge_predict)^2)/length(Y_test)
  ridge_R2= 1 - (sum((Y_test-ridge_predict )^2)/sum((Y_test-mean(Y_test))^2))
  cat('\nRidge R2: ',ridge_R2)
  cat("\nRidge MSE: ",ridge_MSE)
  cat("\n")


  ridge_R2_mean=ridge_R2_mean+ridge_R2
  ridge_MSE_mean=ridge_MSE_mean+ridge_MSE
  lasso_R2_mean=lasso_R2_mean+lasso_R2
  lasso_MSE_mean=lasso_MSE_mean+lasso_MSE
}

cat('\nRidge R2 (mean): ',ridge_R2_mean/10)
cat("\nRidge MSE (mean): ",ridge_MSE_mean/10)
cat('\nLasso R2 (mean): ',lasso_R2_mean/10)
cat("\nLasso MSE (mean): ",lasso_MSE_mean/10)



lasso_lambda = cv.glmnet(y=Y_FULL,x=X_FULL,alpha = 1)
lasso_model = glmnet(X_FULL,Y_FULL,alpha=1,lambda=lasso_lambda$lambda.min)
lasso_coef = which(as.data.frame(as.matrix(coef(lasso_model))) > 0, arr.ind = TRUE)
pred.protein = predict(lasso_model,s=lasso_lambda$lambda.min,newx=data.test)

#lasso_coef shows importance of the individual variables 
X_corr = data.train[,c("x1","x2","x3","x4","x5")]
correlation_selected = cor(na.omit(X_corr))

save(pred.protein, file = "prot_pred.RData")
