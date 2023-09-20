library(MASS)
library(devtools)
library(ggbiplot)
library(e1071)
library(scatterplot3d)
library(ploygon)
library(magrittr)

#fread MC
#data.table
data = read.csv("~/Desktop/carpet/data/data100.csv", header = T, sep = ",")[, -1]
data = as.matrix(data)
data_sd = apply(data, 2, sd)
k = function( x) which(comb_index[1,]==x[1] & comb_index[2,]==x[2])
kinv = function(x) comb_index[, x]

# w_true = diag(c(1.3888, 0.2222, 0.00826))
# C = 1.02
 
w_true = diag(c(2, 2, 1))
C = 10.16

mu = apply(data, 2, mean)
Sigma = cov(data); sd = sqrt(diag(Sigma))
X = data
n = dim(X)[1]; m = 10
comb_index = combn(n, 2)

{
  ## distance data
  Xd = t(apply(comb_index, 2, function(x) X[x[1], ]-X[x[2], ]))

  ## PCA
  X_pc = prcomp(X)
  X_pcs = X %*% X_pc$rotation

  Xd_pc = prcomp(Xd)
  Xd_pcs = Xd %*% Xd_pc$rotation
}

{
  #original data: PCA
  par(mfrow=c(2, 2))
  plot(X_pc)
  plot(X_pcs[,1], X_pcs[,2], xlab="PC1", ylab="PC2", type = "p", cex = 0.3)
  
  #distance data: PCA (not linear transform yet)
  plot(Xd_pc)
  plot(Xd_pcs[,1], Xd_pcs[,2], xlab="PC1", ylab="PC2", type = "p", cex = 0.3)
}

# ## true responses
Y_or = apply(Xd, 1, function(x) t(x) %*% w_true %*% x )
Y = sign(C-Y_or); Y_1 = which(Y==1); Y_2 = which(Y==-1)

# algorithm
rowselection = function(X, m, sam=NULL, M=NULL, plot=T, Y_true=NULL,
                        kernel=NULL, cost=NULL, coef=NULL, degree=NULL, gamma=NULL){
  if(kernel=="radial"){
    dim = 10
    if(is.null(degree)) degree = 3
    if(is.null(gamma)) gamma = 1 / ncol(X)
    if(is.null(coef0)) coef0 = 0
    if(is.null(cost)) cost = 10
    phi=function(x, D){
      vector = 0:(D-1)
      Xs = as.matrix(x) %>% apply(1, function(e) (2*gamma*e^2)^vector / factorial(vector))
      output = apply(Xs, 1, prod) 
      return(output)
    }
    Xs = X %>% apply(1, phi, D=10)
  }else{
    if(is.null(M)) M = diag(dim(X)[2])
    ## linear transformation via eigen-decomposition
    save = eigen(M)
    W = t(save$vectors %*% diag(sqrt(save$values)))
    Xs = X%*%t(W)
  }
  ## Calculate the pairwise Euclidean distance
  ## Xs is the new feature space (standardized, or hyper feature space)
  dist_matrix = as.matrix(dist(Xs, diag = T, method = "euclidean"))
  
  # greedy algorithm for minimum full-connected distance
  greedy_forward = function(sam, Xs, dist_matrix2, m){
    for( j in length(sam):(m-1)){
      dist_matrix2[sam, sam] = Inf
      if(j==1){ 
        add = which.min(dist_matrix2[sam, ]) 
      }else{
        add = order(apply(dist_matrix2[sam, ], 2, sum))[1]
      }
      sam = c(sam, add)
    }
    d = dist(Xs[sam,])
    dis1 = sum(d)
    #dis2 = sum(d^2)
    output = c(sort(sam), dis1)
    return(output)
  }
  
  if(is.null(sam)){
    rec = t(sapply(1:n, greedy_forward, m=m, Xs=Xs, dist_matrix2=dist_matrix))
    sam_min= which.min(rec[,m+1])
    temp = rec[sam_min ,m+1]
    sam = rec[sam_min, -(m+1)]
    sam = sort(sam)
    names(sam)=NULL
  }else{
    sam = greedy_forward(sam=sam, m=m, Xs=Xs, dist_matrix2=dist_matrix)[-(m+1)]  
  }
  # result plot
  if(plot == T)
  {
    par(mfrow=c(1, 2))
    plot( X_pcs[, 1:2], col=1, cex = 0.5, main="scatter plot of S" )
    text(X_pcs[sam, 1]+1, X_pcs[sam, 2], sam, cex = 0.3)
    points( X_pcs[sam, 1:2], col=2, pch=16 )
    index = apply( combn(sam, 2), 2, k)
    
    plot(Xd_pcs[,1:2], col=1, cex=0.1, xlim=c(max(Xd_pcs[index,1]), min(Xd_pcs[index,1]))*1.5,
         ylim=c(max(Xd_pcs[index,2]), min(Xd_pcs[index,2]))*1.5, main="scatter plot and C(S,2)")
    if(is.null(Y_true)){
      points( Xd_pcs[index, 1:2], col=2, pch=16, cex = 0.6)
    }else{
      match = index[Y_true[index]==1]
      points( Xd_pcs[match, 1:2], col=2, pch=16, cex = 0.7)
      unmatch =  index[Y_true[index]==-1]
      points( Xd_pcs[unmatch, 1:2], col=4, pch=16, cex = 0.7)
    }
    par(mfrow=c(1, 1))
  }
  output = list(sam=sam, Xs=Xs, phi=phi, )
  return(sam)
}
rulelearning=function(X, sam, w_true, C, Y_true, plot=T, M=NULL, phi=NULL, 
                      kernel="linear", scale = F, cost=NULL, coef=NULL, degree=NULL, gamma=NULL){
  m = length(sam)
  
  # responses
  index = apply( combn(sam, 2), 2, k)
  Y_des = Y_true[index]; Y_1_des = which(Y_des==1); Y_2_des = which(Y_des==-1)
  # Y_des =Y_true[index] = = sign(C-Y_or_des)
  
  # original X
  X_des = X[sam, ]
  
  # feature difference
  Xd_des = t(apply(combn(m, 2), 2, function(x) (X_des[x[1],] - X_des[x[2],])))
  
  # transfer to feature space
  if(!is.null(M)){ 
    ## linear transformation via eigen-decomposition
    save = eigen(M)
    W = t(save$vectors %*% diag(sqrt(save$values)))
    phi = function(x) x%*%t(W)
  }else if(is.null(phi)){
    phi = function(x) x
  }
  dat_des = data.frame( y=as.factor(Y_des), x=phi(Xd_des) )
  
  
  {
    if(is.null(degree)) degree = 3
    if(is.null(gamma)) gamma = 1 / ncol(X)
    if(is.null(coef0)) coef0 = 0
    if(is.null(cost)) cost = 10
    
    X_svm_des = svm(y~ ., data=dat_des, scale=F, kernel = kernel, cost=cost, gamma=gamma, coef0=coef0)
    ## watch out phi !
    predd = (as.numeric(predict(X_svm_des, data.frame(x=phi(Xd)) ))*2-3)
    table(predd, Y_true)
  }
  
  # summary(X_svm_des)
  beta = t(X_svm_des$coefs) %*% (X_svm_des$SV)
  beta0 = X_svm_des$rho
  a = cbind(beta/beta0, 1)
  
  if( plot==T )
  {
    Y_1 = which(Y==1); Y_2 = which(Y==-1)
    P_1 = which(predd==1); P_2 = which(predd==-1)
    
    # index of choosen sample
    index = apply( combn(sam, 2), 2, k)
    plot(Xd_pcs[,1:2], col=1, cex=0.1, xlim=c(max(Xd_pcs[index,1]), min(Xd_pcs[index,1]))*3,
         ylim=c(max(Xd_pcs[index,2]), min(Xd_pcs[index,2]))*3, main="scatter plot and C(S,2)")
    
    YP11 = intersect(Y_1, P_1)
    YP12 = intersect(Y_1, P_2)
    YP21 = intersect(Y_2, P_1)
    YP22 = intersect(Y_2, P_2)
    
    points( Xd_pcs[YP11, 1:2], col=5, pch=16, cex = 0.5)
    points( Xd_pcs[YP12, 1:2], col=6, pch=16, cex = 0.7)
    points( Xd_pcs[YP21, 1:2], col=7, pch=16, cex = 0.7)
    points( Xd_pcs[YP22, 1:2], col=8, pch=16, cex = 0.7)
    
    match = index[Y_true[index]==1]
    points( Xd_pcs[match, 1:2], col=2, pch=16, cex = 0.7)
    unmatch =  index[Y_true[index]==-1]
    points( Xd_pcs[unmatch, 1:2], col=4, pch=16, cex = 0.7)
  }
  ## after linear transformation, the desicion rule would be different
  ## Gonna try
  output = list( table(predd, Y_true), rbind(a*C,c(diag(w_true), C)) )
  return(output)
}

# full data - result of SVM
result_full1 = rulelearning(X=X, sam=1:n, w_true = w_true, C=C, Y_true=Y, plot=T, 
                            kernel = "linear", phi = function(x) x^2)
result_full2 = rulelearning(X=X, sam=1:n, w_true = w_true, C=C, Y_true=Y, plot=T, 
                            kernel = "radial")

# row selection
m=10
M_maha = solve(cov(X))
M_std = diag(1/apply(X, 2, var))
sam_maha = rowselection(X=X, m=m, M=M_maha, Y_true=Y)
sam_std = rowselection(X=X, m=m, M=M_std, Y_true=Y)
sam_mix = rowselection(X=X, m=m, M=(M_maha+M_std)/2, , Y_true=Y)

print(rbind(sam_maha, sam_std, sam_mix))

# algorithm 2
# sam = intersect(sam_maha, sam_std)
# sam_maha = rowselection(X, m=m, sam=sam, M=M_maha)
# sam_std = rowselection(X, m=m, sam=sam, M=M_std)

# bagging 
# adaboost: robust, lack of weakness


result_maha = rulelearning(X=X, sam=sam_maha, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T)
result_std = rulelearning(X=X, sam=sam_std, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T)
result_mix = rulelearning(X=X, sam=sam_mix, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T)

# learning decision rule
# try X=Xs
result_maha = rulelearning(X=X, sam=sam_std, M=M_maha, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T, cost = 1e2)
result_std = rulelearning(X=X, sam=sam_std, M=M_std, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T)
result_mix = rulelearning(X=X, sam=sam_std, M=M_mix, w_true = w_true, C=C, Y_true=Y, kernel="radial", plot=T)
