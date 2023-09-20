library(MASS)
library(devtools)
library(ggbiplot)
library(e1071)
library(scatterplot3d)

data = read.csv("~/Desktop/carpet/data/data100.csv", header = T, sep = ",")[, -1]
data = as.matrix(data)
data_sd = apply(data, 2, sd)
k = function( x) which(comb_index[1,]==x[1] & comb_index[2,]==x[2])
kinv = function(x) comb_index[, x]

# set.seed(2017)
# w_true = c(1.3888, 0.2222, 0.00826)
# C = 1.02

w_true = c(2, 2, 1)
C = 10.16


mu = apply(data, 2, mean)
Sigma = cov(data); sd = sqrt(diag(Sigma))
X = data
n = dim(X)[1]; m = 10
comb_index = combn(n, 2)

{
  ## distance data
  Xd = t(apply(comb_index, 2, function(x) (X[x[1], ]-X[x[2], ])^2))
  
  ## PCA
  X_pc = prcomp(X)
  X_pcs = X %*% X_pc$rotation
  
  Xd_pc = prcomp(Xd)
  Xd_pcs = Xd %*% Xd_pc$rotation
}

#original data: PCA
par(mfrow=c(1, 2))
plot(X_pc)
plot(X_pcs[,1], X_pcs[,2], xlab="PC1", ylab="PC2", type = "p", cex = 0.3)


#distance data: PCA
plot(Xd_pc)
plot(Xd_pcs[,1], Xd_pcs[,2], xlab="PC1", ylab="PC2", type = "p", cex = 0.3)



Y_or = (Xd^2 %*% (w_true))
Y = sign(C-Y_or); Y_1 = which(Y==1); Y_2 = which(Y==-1)
length(Y_1)
dat = data.frame(y=as.factor(Y), x=Xd)

X_svm = svm(y~ ., data=dat, kernel = "linear", cost=1e10, scale=F)
summary(X_svm)
SV = X_svm$index
select = as.numeric(comb_index[,SV])
{
  par(mfrow=c(1, 2))
  plot(X_pcs[,1:2], cex = 0.5, xlab = "PC1", ylab="PC2")
  text(X_pcs[select,1:2], X_pcs[select, 2], select, col='2')
  
  
  plot( Xd_pcs[, 1:2], col=1, cex=0.2, xlim=c(0, 6), ylim = c(-6, 1.2), xlab = "X1", ylab="X2")
  points( Xd_pcs[Y_1, 1], Xd_pcs[Y_1, 2], col=2, cex=0.2 )
  
  beta = t(X_svm$coefs) %*% (X_svm$SV)
  beta0 = X_svm$rho
  a = cbind(beta/beta0, 1)
  print(a*C)
  
  select_dist = t(apply(combn(select, 2), 2, function(x) (X[x[1], ]-X[x[2], ])^2 ))
  points( select_dist %*% Xd_pc$rotation[,1:2], col=4, cex = .8, pch=16)
  points( X_svm$SV %*% Xd_pc$rotation[,1:2], col=3, cex=1, pch=16)
  
  par(mfrow=c(1, 1))
}

## algorithm
#algorithm
m=10
## Mahalanobis
M = solve(cov(X))
save = eigen(M)
W = t(save$vectors %*% diag(sqrt(save$values)))
W = diag(1/apply(X, 2, sd))

Xs = X%*%t(W)

dist_matrix = as.matrix(dist(Xs, diag = T, method = "euclidean"))
rec = matrix(0, n, m+1)
for(i in 1:n){
  sam = i
  dist_matrix2=dist_matrix
  for( j in 1:(m-1)){
    dist_matrix2[sam, sam] = Inf
    if(j==1) add = which.min(dist_matrix2[sam, ]) else{
      add = order(apply(dist_matrix2[sam, ], 2, sum))[1]
    }
    sam = c(sam, add)
    sam
  }
  d = dist(Xs[sam,])
  dis1 = sum(d)
  dis2 = sum(d^2)
  rec[i, ] = c(sam, dis1)
}
sam_min= which.min(rec[,m+1])
temp = rec[,m+1]
sam = rec[sam_min, -(m+1)]
sam = sort(sam)
print(sam)

par(mfrow=c(1, 2))
plot( X_pcs[, 1:2], col=1, cex = 0.5, main="scatter plot and S" )
#text(X_pcs[sam, 1]+1, X_pcs[sam, 2], sam, cex = 0.8)
points( X_pcs[sam, 1:2], col=2, pch=16 )
plot(Xd_pcs[,1:2], col=1, cex=0.1, xlim=c(0, 200), ylim=c(-10, 5), main="scatter plot and C(S,2)")
#index = apply( combn(sam, 2), 2, k)
points( Xd_pcs[index, 1:2], col=2, pch=16, cex = 0.5)


## learning decision rule
X_des = X[sam,]
Xd_des = t(apply(combn(m, 2), 2, function(x) (X_des[x[1],] - X_des[x[2],])^2 ))

#t(apply(combn(m, 2), 2, function(x) (X_des[x[1], ]-X_des[x[2], ])^2 )) == Xd_des


Y_or_des = (Xd_des %*% (w_true))
Y_des = sign(C-Y_or_des); Y_1_des = which(Y_des==1); Y_2_des = which(Y_des==-1)

#index = apply(combn(sam, 2), 2,k)
#Y_des * Y[index]
table(Y_des)
dat_des = data.frame(y=as.factor(Y_des), x=Xd_des)

X_svm_des = svm(y~ ., data=dat_des, kernel = "linear", scale=F, cost=10000)
summary(X_svm_des)
{
  plot( Xd_pcs[, 1:2], col=1, cex=0.2, xlim=c(0, 100), ylim = c(-20, 1.2), xlab = "X1", ylab="X2")
  points( Xd_pcs[Y_1, 1:2], col=2, cex=0.2)
  
  points( Xd_des[Y_1_des, ] %*% Xd_pc$rotation, col=3, cex=1, pch=16 )
  points( Xd_des[Y_2_des, ] %*% Xd_pc$rotation, col=4, cex=1, pch=16 )
  points( X_svm_des$SV %*% Xd_pc$rotation, col=5, cex=1, pch=16 )
}


beta = t(X_svm_des$coefs) %*% (X_svm_des$SV)
beta0 = X_svm_des$rho
a = cbind(beta/beta0, 1)
print( rbind(a*C,c(w_true, C)) )

predd = (as.numeric(predict(X_svm_des, data.frame(x=Xd) ))*2-3)
wrong = which(predd*Y== -1)
table(predd, Y)
points( Xd_pcs[wrong, 1:2], col=6, cex=1, pch=16 )

