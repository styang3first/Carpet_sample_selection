library(MASS)
library(devtools)
library(ggbiplot)
library(e1071)
library(scatterplot3d)
library(ploygon)
library(magrittr)
library(randtoolbox)
repeat_permut=function(j, p){
  save = rbind(0, combn(j+p-1, p-1))
  save = apply(save, 2, diff)-1
  if(p==2) save = t(save)
  save = rbind(save, j-colSums(save))
  return(save)
}
expand=function(x, j){
  p = length(x)
  count = repeat_permut(j=j, p=p)
  # column product
  output = apply(count, 2, function(v) x^v/sqrt(factorial(v))) %>% apply(2, prod)
  return(output)
}
phi_kernel = function(x, gamma, D){
  v = 0:D
  p = length(x)
  feature_space = apply(as.matrix(v), 1, expand, x=x*sqrt(2*gamma))
  output = NULL
  for(i in (v+1)) output = c(output, feature_space[[i]])
  output = output * exp(-gamma*sum(x^2))
  return(output)
}
greedy_forward = function(sam, Xs, dist_matrix2, m){
  copy = dist_matrix2
  for( j in length(sam):(m-1)){
    dist_matrix2[sam, sam] = Inf
    if(j==1){ 
      add = which.min(dist_matrix2[sam, ]) 
    }else{
      add = order(apply(dist_matrix2[sam, ], 2, sum))[1]
    }
    sam = c(sam, add)
  }
  d = as.dist(copy[sam, sam])
  dis1 = sum(d)
  output = c(sort(sam), dis1)
  return(output)
}
visualplot = function(x, index, Y=NULL, lim=NULL, pcaplot=F, main){
  x_pc = prcomp(x)
  if(pcaplot){ 
    par(mfrow=c(1, 2))
    plot(x_pc)
  }
  x_pcs = x %*% x_pc$rotation
  
  if(!is.null(lim)){ 
    xlim=range(x_pcs[index,1])*1.5
    ylim=range(x_pcs[index,2])*1.5
  }else{
    xlim=ylim=NULL
  }
  
  plot(x_pcs[,1:2], col=1, cex=0.1, xlim=xlim, ylim=ylim, main=main)
  if(is.null(Y)){
    points( x_pcs[index, 1:2], col=2, pch=16, cex = 0.6)
  }else{
    match = index[Y[index]==1]
    points( x_pcs[match, 1:2], col=2, pch=16, cex = 0.7)
    unmatch =  index[Y[index]==-1]
    points( x_pcs[unmatch, 1:2], col=4, pch=16, cex = 0.7)
  }
}
rowselection = function(X, m, sam=NULL, M=NULL, plot=T, kernelplot=F, Y=NULL, distance="linear",
                        kernel="radial", cost=NULL, coef0=NULL, degree=NULL, gamma=NULL){
  # distance data
  X_pc = prcomp(X)
  X_pcs = X %*% X_pc$rotation
  
  if(is.null(M)) M = diag(ncol(X))
  if(is.null(degree)) degree = 3
  if(is.null(gamma)) gamma = 1 / ncol(X)
  if(is.null(coef0)) coef0 = 0
  if(is.null(cost)) cost = 10
  
  ## linear transformation via eigen-decomposition
  save = eigen(M)
  W = t(save$vectors %*% diag(sqrt(save$values)))
  
  ## Xs is the standardized feature space
  Xs = X%*%t(W)
  
  ## Calculate the pairwise distance, Euclidiea, RBF
  if(distance=="linear"){
    dist_matrix = as.matrix(dist(Xs, diag = T, method = "euclidean"))
    phi_kernel = function(x) x
  }else if(distance=="radial"){
    dist_matrix = as.matrix(dist(Xs, diag = T, method = "euclidean"))
    dist_matrix = 2-exp(-gamma*dist_matrix)^2
  }
  
  ## greedy algorithm for minimum full-connected distance
  if(is.null(sam)){
    rec = t(sapply(1:n, greedy_forward, m=m, dist_matrix2=dist_matrix))
    sam_min= which.min(rec[,m+1])
    temp = rec[sam_min ,m+1]
    sam = rec[sam_min, -(m+1)]
    sam = sort(sam)
    names(sam)=NULL
    index = apply( combn(sam, 2), 2, k)
  }else{
    sam = greedy_forward(sam=sam, m=m, dist_matrix2=dist_matrix)[-(m+1)]  
    index = apply( combn(sam, 2), 2, k)
  }
  
  ## result plot
  if(plot == T)
  {
    par(mfrow=c(2, 2))
    Xd = t(apply(comb_index, 2, function(x) X[x[1], ]-X[x[2], ]))
    visualplot(x=X, index=sam, main="scatter plot of original data")
    visualplot(x=Xs, index=sam, main="scatter plot of rotated data")
    visualplot(x=Xd%*%t(W), index=index, Y=Y, lim=1, main="scatter plot of distance data")
    
    # this is not necessary, since the num of dimension is infinity
    if(kernelplot==T){
      Xd_fs = apply(Xd%*%t(W), 1, phi_kernel, gamma=gamma, D=10) %>% t
      visualplot(x=Xd_fs, index=index, Y=Y, lim=1, main="distance data in the feature space")
    }
    par(mfrow=c(1, 1))
  }
  output = list(Y=Y, X=X, Xd=Xd, W=W, sam=sam, index=index, kernel=kernel, 
                cost=cost, coef0=coef0, degree=degree, gamma=gamma)
  # Y: responses
  # X: original data
  # Xd: distance data w/o rotation
  
  # W: roration matrix
  # Xs: manually rotated data
  # sam: selected subsample
  # index: the corresponding index of sam in Xd
  # kernel: for SVM
  return(output)
}
selection_object = rowselection(X=X, m=m, M=M_std, Y=Y, distance="radial", gamma=g)

rulelearning=function(selection_object, phi_manual=1, plot=T,
                      scale=F, cost=NULL, coef0=NULL, degree=NULL, gamma=NULL, cross=0){
  attach(selection_object)
  if(is.null(cost)) cost = selection_object$cost
  if(is.null(coef0)) coef0 = selection_object$coef0
  if(is.null(degree)) degree = selection_object$degree
  if(is.null(gamma)) gamma = selection_object$gamma
  m = length(sam)
  
  if(phi_manual=="rotate"){ 
    phi_manual = function(x) x%*%t(W)
  }else{
    phi_manual = function(x) x
  }
  Xd_des = Xd[index,]
  Y_des = Y[index]
  dat_des = data.frame(y=as.factor(Y_des), x=phi_manual(Xd_des))
  
  
  ## conduct cross-validation
  if(cross>0){
    Y1=which(Y_des==1); Y2=which(Y_des==-1)
    l1 = length(Y1); l2=length(Y2)
    cross_index = c(sample(Y1, l1), sample(Y2, l2)) %>% matrix(ncol=cross, byrow=T)
    
    ## grid_search
    C.try=2^(-10:10)
    gamma.try = 2^(0:14)
    parameter = expand.grid(C=C.try, g=gamma.try)
    parameter = rbind(parameter, c(10, 1/3))
    error_compute = function(test_index, type1, data, scale, kernel, cost, gamma, coef0){
        type2 = 1-type1
        training = data[-test_index,]
        testing = data[test_index,]
        svm_object = svm(y~ ., data=training, scale=scale, kernel = kernel, cost=cost, gamma=gamma, coef0=coef0)
        pred = predict(svm_object, testing) %>% as.numeric * 2-3
        Y = testing$y %>% as.numeric * 2-3
        P1 = which(pred==1); P2 = which(pred==-1)
        Y1 = which(Y==1); Y2 = which(Y==-1)
        
        YP11 = intersect(Y1, P1) %>% length
        YP12 = intersect(Y1, P2) %>% length
        YP21 = intersect(Y2, P1) %>% length
        YP22 = intersect(Y2, P2) %>% length
        
        error = type1*YP21/(YP11+YP21)# + type2*YP12/(YP12+YP22)
        return(error)
      }
    cross_validation = function(data, type1=1, cross_index, scale=NULL, kernel=NULL, cost=NULL, gamma=NULL, coef0=NULL){
      loss = apply(cross_index, 2, function(x){ 
        error_compute(test_index=x, type1=1, data=data, 
                      scale=scale, kernel=kernel, cost=cost, gamma=gamma, coef0=coef0)
      })
      return(sum(loss))
    }
    
    para_tune = apply(parameter, 1, function(x) cross_validation(data=dat_des, type1=1, 
                        cross_index=cross_index, scale=T, kernel=kernel,
                        cost=x[1], gamma=x[2], coef0=coef0)
    )
    cv = which.min(para_tune)  
    cost = parameter[cv, 1]
    gamma = parameter[cv, 2]
  }
  
  
  {
    X_svm_des = svm(y~ ., data=dat_des, scale=scale, kernel = kernel, cost=cost, gamma=gamma, coef0=coef0)
    ## watch out phi_manual !
    predd = predict(X_svm_des, data.frame(x=phi_manual(Xd)) ) %>% as.numeric *2-3
    # print(table(predd, Y))
    
    dat = data.frame(y=as.factor(Y), x=phi_manual(Xd))
  }
  
  
  
  if( plot==T )
  {
    P1 = which(predd==1); P2 = which(predd==-1)
    Y1 = which(Y==1); Y2 = which(Y==-1)
    
    # index of choosen sample
    YP11 = intersect(Y1, P1)
    YP12 = intersect(Y1, P2)
    YP21 = intersect(Y2, P1)
    YP22 = intersect(Y2, P2)
    
    plot(Xd_pcs[,1:2], col=1, cex=0.1, xlim=range(Xd_pcs[index,1])*2,
         ylim=range(Xd_pcs[index,2])*2, main="scatter plot and C(S,2)")
    
    points( Xd_pcs[YP11, 1:2], col=2, pch=16, cex = 0.5)
    points( Xd_pcs[YP12, 1:2], col=3, pch=16, cex = 0.7)
    points( Xd_pcs[YP21, 1:2], col=4, pch=16, cex = 0.7)
    legend("topleft", c("(Y,Y')=(1, 1)", "(Y,Y')=(1, 2)", "(Y,Y')=(2, 1)", "(Y,Y')=(2, 2)"), pch = rep(16, 4), col=c(2:4, 1))
    
    if(ncol(X)==2) {
      index2 = union( sample(1:4950, 1000), index)
      plot(X_svm_des, data=dat[index2, ], grid=500, xlim=range(Xd[index,1])*2,
           ylim=range(Xd[index,2])*2)
    }
  }
  
  ## after linear transformation, the desicion rule would be different
  ## Gonna try
  output = list( table(predd, Y))
  return(output)
}
#fread MC
#data.table
data = read.csv("~/Desktop/carpet/data/data100.csv", header = T, sep = ",")[, -1]
data = as.matrix(data)
data_sd = apply(data, 2, sd)
k = function(x) which(comb_index[1,]==x[1] & comb_index[2,]==x[2])
kinv = function(x) comb_index[, x]

# w_true = diag(c(1.3888, 0.2222, 0.00826))
# C = 1.02

w_true = diag(c(2, 2, 1))
C = 10.16

mu = apply(data, 2, mean)
Sigma = cov(data); sd = sqrt(diag(Sigma))
X = data
n = dim(X)[1]
comb_index = combn(n, 2)

#original data: PCA
visualplot(X, index=1:100, main="scatter plot", pcaplot=T)
#distance data: PCA (not linear transform yet)
visualplot(Xd, index=1:nrow(Xd), main="scatter plot", pcaplot=T)

# ## true responses
Y_or = apply(Xd, 1, function(x) t(x) %*% w_true %*% x )
Y = sign(C-Y_or); Y_1 = which(Y==1); Y_2 = which(Y==-1)
# row selection algorithm

m=10
M_maha = solve(cov(X))
M_std = diag(1/apply(X, 2, var))
M_mix = (M_maha+M_std)/2
g=1/3
select_maha_linear = rowselection(X=X, m=m, M=M_maha, Y=Y, gamma=g)
select_maha_radial = rowselection(X=X, m=m, M=M_maha, Y=Y, distance="radial", gamma=g)
select_std_linear = rowselection(X=X, m=m, M=M_std, Y=Y, gamma=g)
select_std_radial = rowselection(X=X, m=m, M=M_std, Y=Y, distance="radial", gamma=g)
select_mix_linear = rowselection(X=X, m=m, M=M_mix, Y=Y, gamma=g)
select_mix_radial = rowselection(X=X, m=m, M=M_mix, Y=Y, distance="radial", gamma=g)
select_linear = rowselection(X=X, m=m, Y=Y, gamma=g)
select_radial = rowselection(X=X, m=m, Y=Y, distance="radial", gamma=g)



print(rbind(select_maha_linear$sam, 
            select_maha_radial$sam, 
            select_std_linear$sam, 
            select_std_radial$sam, 
            select_mix_linear$sam, 
            select_mix_radial$sam,
            select_linear$sam, 
            select_radial$sam))

# bagging 
# adaboost: robust, lack of weakness
par(mfrow=c(4, 2))
phi = function(x) x^2
phi = function(x) x
result_maha_linear = rulelearning(select_maha_linear)
result_maha_radial = rulelearning(select_maha_radial)
result_std_linear = rulelearning(select_std_linear)
result_std_radial = rulelearning(select_std_radial)
result_std_linear_scale = rulelearning(select_std_linear, scale=T)
result_std_radial_scale = rulelearning(select_std_radial, scale=T)

result_mix_linear = rulelearning(select_mix_linear)
result_mix_radial = rulelearning(select_mix_radial)
result_linear = rulelearning(select_linear)
result_radial = rulelearning(select_radial)
print(result_maha_linear)
print(result_maha_radial)
print(result_std_linear)
print(result_std_radial)
print(result_mix_linear)
print(result_mix_radial)
print(result_linear)
print(result_radial)



# why svn perform worse after we standardize the distance matrix
result_maha_linear = rulelearning(select_maha_linear, phi_manual="rotate")
result_maha_radial = rulelearning(select_maha_radial, phi_manual="rotate")
result_std_linear = rulelearning(select_std_linear, phi_manual="rotate")
result_std_radial = rulelearning(select_std_radial, phi_manual="rotate")

rulelearning(select_std_radial, phi_manual=1, cost=10000, gamma=1/1.8)

result_mix_linear = rulelearning(select_mix_linear, phi_manual="rotate")
result_mix_radial = rulelearning(select_mix_radial, phi_manual="rotate")
result_linear = rulelearning(select_linear, phi_manual="rotate")
result_radial = rulelearning(select_radial, phi_manual="rotate")

print(result_maha_linear)
print(result_maha_radial)
print(result_std_linear)
print(result_std_radial)
print(result_mix_linear)
print(result_mix_radial)
print(result_linear)
print(result_radial)


## summary:
In row selection:
1. need to rotate X to Xs, or the distance would be unfair
2. Distance is measured by Euclidean or RBF
3. Mahalanobis standarization performed worse. However, the reasone is 
   probably because the true model doesnt have interaction
4. the Bayesian result looks good, but #YP=(-1, 1) is too large
5. Standardize perfoerms good

In rule decision:
1. I dont know why the result is worse if I use phi(x) = x%*%t(W) as predictor
2. use the original Xd is enough
3. probably that is I dont know how to tune SVM?

How to tune?
1. scale data = [-1, 1] or [0, 1]
2. (C, gamma): exponentially grid-search