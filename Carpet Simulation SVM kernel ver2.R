library(MASS)
library(devtools)
library(ggbiplot)
library(e1071)
library(scatterplot3d)
library(ploygon)
library(magrittr)
library(randtoolbox)
library(data.table)

k = function(x, comb_index) which(comb_index[1,]==x[1] & comb_index[2,]==x[2])
k2 = function(x, n, comb_index) -n + (n-0.5)*x[1, ]+x[2, ]-0.5*x[1, ]^2
distdata = function(x, comb_index){
  output = apply(comb_index, 2, function(ind) x[ind[1], ]-x[ind[2], ]) %>% t
  return(output)
}
greedy_forward1 = function(sam, m, n, dist_matrix1){
  while(length(sam)<m){
    rest = (1:n)[-sam]
    if(length(sam)==1){
      cost = dist_matrix1[sam, rest]
    }else{
      cost = colSums(dist_matrix1[sam, rest])
    }
    add = rest[cost%>%which.min]
    sam = c(sam, add)
  }
  sam = sort(sam)
  d = as.dist(dist_matrix1[sam, sam])
  dis1 = sum(d)
  output = c(sam, dis1)
  return(output)
}
greedy_forward2 = function(sam, m, n, dist_matrix2, comb_index){
  print(sam)
  if(length(sam)==1) return("no result")
  while(length(sam)<m){
    rest = (1:n)[-sam]
    cost = sapply(rest, function(a){
      new = c(sam, a) %>% sort
      # sam1_new = apply( combn(new, 2), 2, k, comb_index=comb_index)
      sam1_new = combn(new, 2) %>% k2( n=n, comb_index=comb_index)
      #as.dist(dist_matrix2[sam1_new, sam1_new])
      cost = sum(dist_matrix2[sam1_new, sam1_new])
    })
    add = rest[cost%>%which.min]
    sam = c(sam, add)
  }
  sam = sort(sam)
  sam1 = combn(sam, 2) %>% k2( n=n, comb_index=comb_index)
  dis1 = sum(dist_matrix2[sam1, sam1])
  output = c(sam, dis1)
  return(output)
}


visualplot = function(x, sam, Y=NULL, lim=NULL, pcaplot=F, main){
  x_pc = prcomp(x)
  if(pcaplot){ 
    par(mfrow=c(1, 2))
    plot(x_pc)
  }
  x_pcs = x %*% x_pc$rotation
  
  if(!is.null(lim)){ 
    xlim=range(x_pcs[sam,1])*1.5
    ylim=range(x_pcs[sam,2])*1.5
  }else{
    xlim=ylim=NULL
  }
  
  plot(x_pcs[,1:2], col=1, cex=0.1, xlim=xlim, ylim=ylim, main=main)
  if(is.null(Y)){
    points( x_pcs[sam, 1:2], col=2, pch=16, cex = 0.6)
  }else{
    match = sam[Y[sam]==1]
    points( x_pcs[match, 1:2], col=2, pch=16, cex = 0.7)
    unmatch =  sam[Y[sam]==-1]
    points( x_pcs[unmatch, 1:2], col=4, pch=16, cex = 0.7)
  }
}
# M: data rotating parameter
# plot: selection result plot
# kernelplot: selection result plot in feature space
{
  data = read.csv("~/Desktop/carpet/data/data100.csv", header = T, sep = ",")[, -1]
  data = as.matrix(data)
  data_sd = apply(data, 2, sd)
  w_true = diag(c(2, 2, 1))
  C = 10.16
  mu = apply(data, 2, mean)
  Sigma = cov(data); sd = sqrt(diag(Sigma))
  X = data
  Xd = distdata(X, comb_index=combn(100, 2) )
  Y_or = apply(Xd, 1, function(x) t(x) %*% w_true %*% x )
  Y = sign(C-Y_or); Y_1 = which(Y==1); Y_2 = which(Y==-1)
  m=10
  M = diag(1/apply(X, 2, var))
}
rowselection = function(X, m, sam=NULL, M=NULL, plot=T, Y=NULL, criterion=2,
                        kernel="radial", cost=NULL, coef0=NULL, degree=NULL, gamma=NULL){
  # scaling parameter default
  if(is.null(M)) M = diag(ncol(X))
  
  # kernel parameter default
  if(is.null(degree)) degree = 3
  if(is.null(gamma)) gamma = 1 / ncol(X)
  if(is.null(coef0)) coef0 = 0
  if(is.null(cost)) cost = 10
  
  ## linear transformation via eigen-decomposition
  save = eigen(M)
  W = t(save$vectors %*% diag(sqrt(save$values)))
  
  ## Rotated data
  n = nrow(X)
  Xs = X%*%t(W)
  
  comb_index = combn(n, 2)
  ## Calculate the pairwise distance, Euclidiea, RBF

  dist_matrix1 = as.matrix(dist(Xs, diag = T, method = "euclidean"))
  if(criterion==2){
    Xds = distdata(Xs, comb_index=comb_index)
    dist_matrix2 = as.matrix(dist(Xds, diag = T, method = "euclidean"))
  }
  
  if(kernel=="radial"){
    dist_matrix1 = 2-exp(-gamma*dist_matrix1)^2
    if(criterion==2){
      dist_matrix2 = 2-2*exp(-gamma*dist_matrix2)^2
      }
  }
  
  ## greedy algorithm for minimum full-connected distance
  # step 1
  if(is.null(sam)){
    rec = t(sapply(1:n, greedy_forward1, m=m, n=n, dist_matrix1=dist_matrix1))
    sam_min= which.min(rec[,m+1])
    temp = rec[sam_min ,m+1]
    sam = rec[sam_min, -(m+1)]
    sam = sort(sam)
    sam1 = combn(sam, 2) %>% k2(n=n, comb_index=comb_index)
    # any case: 10 16 34 38 61 67 73 75 78 93
  }else{
    sam = greedy_forward1(sam=sam, m=m, n=n, dist_matrix1=dist_matrix1)[-(m+1)]  
    sam1 = combn(sam, 2) %>% k2(n=n, comb_index=comb_index)
  }
  # step 2
  if(criterion==2){
    rec = t(apply(combn(sam, 2), 2, greedy_forward2, m=m, n=n, dist_matrix2=dist_matrix2, comb_index=comb_index))
    # greedy_forward2(sam=1:2, m=m, n=n, dist_matrix2=dist_matrix2, comb_index=comb_index)
    sam_min= which.min(rec[,m+1])
    temp = rec[sam_min ,m+1]
    sam = rec[sam_min, -(m+1)]
    sam = sort(sam)
    sam1 = combn(sam, 2) %>% k2(n=n, comb_index=comb_index)
    # linear:            10  16  34  38  46  61  73  78  85 100
    # radial(gamma=5):   10  16  34  38  61  67  73  78  85 100
    # radial(gamma=1/3): 10  16  34  38  46  61  73  78  85 100
  }
  ## result plot
  if(plot == T)
  {
    par(mfrow=c(1, 2))
    Xd=distdata(X, comb_index=comb_index )
    
    visualplot(x=X, sam=sam, main="scatter plot of original data")
    #visualplot(x=Xs, sam=sam, main="scatter plot of rotated data")
    visualplot(x=Xd, sam=sam1, Y=Y, lim=1, main="scatter plot of distance data")
    
    par(mfrow=c(1, 1))
  }
  output = list(Y=Y, X=X, W=W, sam=sam, sam1=sam1, comb_index=comb_index, kernel=kernel, 
                cost=cost, coef0=coef0, degree=degree, gamma=gamma)
  # Y: responses
  # X: original data
  # W: roration matrix
  # Xs: rotated data
  # sam: selected subsample
  # sam1: the corresponding index of sam in Xd
  # sam2: the corredponding index of sam in Xd*
  # kernel: used kernel
  return(output)
}
rulelearning=function(selection_object, plot=T,
                      kernel=NULL, scale=F, cost=NULL, coef0=NULL, degree=NULL, gamma=NULL, cross=0){
  attach(selection_object)
  if(is.null(kernel)) kernel = selection_object$kernel
  if(is.null(cost)) cost = selection_object$cost
  if(is.null(coef0)) coef0 = selection_object$coef0
  if(is.null(degree)) degree = selection_object$degree
  if(is.null(gamma)) gamma = selection_object$gamma
  m = length(sam)
  Xd = distdata(X, comb_index)
  Xd_des = Xd[sam1,]
  Y_des = Y[sam1]
  dat_des = data.frame(y=as.factor(Y_des), x=Xd_des)
 
  # conduct cross-validation
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
    predd = predict(X_svm_des, data.frame(x=Xd) ) %>% as.numeric *2-3
    # print(table(predd, Y))
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
    
    
    Xd_pc = prcomp(Xd)
    Xd_pcs = Xd %*% Xd_pc$rotation
    plot(Xd_pcs[,1:2], col=1, cex=0.1, xlim=range(Xd_pcs[sam1,1])*2,
         ylim=range(Xd_pcs[sam1,2])*2, main="scatter plot and C(S,2)")
    
    points( Xd_pcs[YP11, 1:2], col=2, pch=16, cex = 0.5)
    points( Xd_pcs[YP12, 1:2], col=3, pch=16, cex = 0.7)
    points( Xd_pcs[YP21, 1:2], col=4, pch=16, cex = 0.7)
    legend("topleft", c("(Y,Y')=(1, 1)", "(Y,Y')=(1, 2)", "(Y,Y')=(2, 1)", "(Y,Y')=(2, 2)"), pch = rep(16, 4), col=c(2:4, 1))
  }
  output = list( table(predd, Y) )
  return(output)
}
#fread MC
#data.table
data = read.csv("~/Desktop/carpet/data/data100.csv", header = T, sep = ",")[, -1]
data = as.matrix(data)
data_sd = apply(data, 2, sd)

# w_true = diag(c(1.3888, 0.2222, 0.00826))
# C = 1.02

w_true = diag(c(2, 2, 1))
C = 10.16

mu = apply(data, 2, mean)
Sigma = cov(data); sd = sqrt(diag(Sigma))
X = data

#original data: PCA
visualplot(X, sam=1:nrow(X), main="scatter plot", pcaplot=T)
#distance data: PCA (not linear transform yet)
visualplot(Xd, sam=1:nrow(Xd), main="scatter plot", pcaplot=T)

# ## true responses
comb_index = combn(nrow(X), 2)
Xd = distdata(X, comb_index = comb_index)
Y_or = apply(Xd, 1, function(x) t(x) %*% w_true %*% x )
Y = sign(C-Y_or); Y_1 = which(Y==1); Y_2 = which(Y==-1)
table(Y)
# row selection algorithm

m=10
M_maha = solve(cov(X))
M_std = diag(1/apply(X, 2, var))
M_mix = (M_maha+M_std)/2

select_maha_linear1 = rowselection(X=X, m=m, M=M_maha, Y=Y, kernel="linear", criterion=1)
select_maha_linear2 = rowselection(X=X, m=m, M=M_maha, Y=Y, kernel="linear", criterion=2)
select_std_linear1 = rowselection(X=X, m=m, M=M_std, Y=Y, kernel="linear", criterion=1)
select_std_linear2 = rowselection(X=X, m=m, M=M_std, Y=Y, kernel="linear", criterion=2)


select_maha_radial1 = rowselection(X=X, m=m, M=M_maha, Y=Y, kernel="radial", criterion=1, gamma=3)
select_maha_radial2 = rowselection(X=X, m=m, M=M_maha, Y=Y, kernel="radial", criterion=2, gamma=3)
print(rbind(select_maha_linear1$sam, 
            select_maha_linear2$sam, 
            select_maha_radial1$sam, 
            select_maha_radial2$sam))

rulelearning(select_maha_linear1)
rulelearning(select_maha_linear2)
rulelearning(select_maha_radial1)
rulelearning(select_maha_radial1)



select_std_linear = rowselection(X=X, m=m, M=M_std, Y=Y, gamma=g)
select_std_radial = rowselection(X=X, m=m, M=M_std, Y=Y, distance="radial", gamma=g)
select_mix_linear = rowselection(X=X, m=m, M=M_mix, Y=Y, gamma=g)
select_mix_radial = rowselection(X=X, m=m, M=M_mix, Y=Y, distance="radial", gamma=g)
select_linear = rowselection(X=X, m=m, Y=Y, gamma=g)
select_radial = rowselection(X=X, m=m, Y=Y, distance="radial", gamma=g)
selection_object1 = rowselection(X=X, m=m, M=M, Y=Y, kernel="radial", gamma=1/3)
selection_object2 = rowselection(X=X, m=m, M=M, Y=Y, kernel="radial", gamma=5)
rulelearning(selection_object=selection_object1, gamma=5)
rulelearning(selection_object=selection_object2)



      , 
            select_std_linear$sam, 
            select_std_radial$sam, 
            select_mix_linear$sam, 
            select_mix_radial$index,
            select_linear$index, 
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
result_maha_linear = rulelearning(select_maha_linear)
result_maha_radial = rulelearning(select_maha_radial)
result_std_linear = rulelearning(select_std_linear)
result_std_radial = rulelearning(select_std_radial)

rulelearning(select_std_radial, cost=10000, gamma=1/1.8)

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





ref: https://rpubs.com/skydome20/R-Note14-SVM-SVR

require("mlbench")
data(Glass, package="mlbench")
data = Glass

# tune cost and gamma in SVM(soft-margin)
tune.model = tune(svm,
                  Type~.,
                  data=data,
                  kernel="radial", # RBF kernel function
                  range=list(cost=10^(-1:2), gamma=c(.5,1,2))# 調參數的最主要一行
)
summary(tune.model)
plot(tune.model)
tune.model$best.model