
## CrossMST

getR1R2 = function(E, treated.index){
  R1 = R2 = 0
  for (i in 1:nrow(E)){
    e1 = is.na(match(E[i,1],treated.index))
    e2 = is.na(match(E[i,2],treated.index))
    if ((!e1) && (!e2))  R1 = R1 + 1
    if (e1 && e2)  R2 = R2 + 1
  }
  return(list(R1=R1, R2=R2))
}


## function CrossMST calculates the test statistic R_M and p-value for the new test
## if perm>0, also calculate permutation p-value

CrossMST = function(distM,treated.index,perm=0){
  n = length(treated.index)
  N = dim(distM)[1]
  if (N!=(2*n)){
    return("Control index wrong!")
  }
  case.index = (1:N)[-treated.index]
  reorder = c(treated.index,case.index)
  distM2 = distM[reorder,reorder]
  E = mstree(as.dist(distM2))
  Ebynode = vector("list", N)
  for (i in 1:N) Ebynode[[i]] = rep(0,0)
  for (i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]], E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]], E[i,1])
  }
  nE = nrow(E)
  nodedeg = rep(0,N)
  for (i in 1:N) nodedeg[i] = length(Ebynode[[i]])
  nEi = sum(nodedeg*(nodedeg-1))  # pair of nodes sharing a node * 2
  mu = (N-2)/4
  V = (nEi*N*(N-4)/(N-1) - (N-2)*(N-6))/16/(N-3)
  V12 = (3*(N-2)-nEi*N/(N-1))*(N-2)/16/(N-3)
  temp = getR1R2(E,1:n)
  Z1 = (temp$R1-mu)/sqrt(V)
  Z2 = (temp$R2-mu)/sqrt(V)
  rho = V12/V
  x = max(Z1,Z2)
  # p1 = 2*pnorm(-x)
  p2 =  1-pmvnorm(upper=rep(x,2),mean=rep(0,2),corr=matrix(c(1,rho,rho,1),2))[1]
  if (perm<=0){
    return(list(test.stat.R=max(temp$R1,temp$R2), test.stat.Z=x, pval.appr=p2))
  }else{
    stat = rep(0,perm)
    for (i in 1:perm){
      temp2 = getR1R2(E,sample(1:N,n))
      stat[i] = max(temp2$R1,temp2$R2)
    }
    p3 = length(which(stat>=max(temp$R1,temp$R2)))/perm
    return(list(test.stat.R=max(temp$R1,temp$R2), test.stat.Z=x, pval.appr=p2,  pval.perm=p3))
  }
}



## function CrossNN calculates the test statistic and p-value for the new test
## if perm>0, also calculate permutation p-value

CrossNN = function(distM,treated.index,perm=0){
  n = length(treated.index)
  N = dim(distM)[1]
  if (N!=(2*n)){
    return("Control index wrong!")
  }
  case.index = (1:N)[-treated.index]
  reorder = c(treated.index,case.index)
  distM2 = distM[reorder,reorder]
  diag(distM2) = max(distM)+10
  nearest = apply(distM2, 1, which.min)
  n = length(treated.index)
  D11 = length(which(nearest[1:n]<=n))
  D22 = length(which(nearest[(n+1):(2*n)]>n))
  D12 = n-D11 # number of nodes in sample 1 whose nearest neighbor is from sample 2
  D21 = n-D22
  
  ### pick one arrow, pick another arrow, 2n*(2n-1) possibilities
  
  ## the number of possibilities that the two arrows share the endpoint
  temp = as.vector(table(nearest))
  k = max(temp)
  share = 0
  if (k>1){
    for (i in 2:k){
      share = share + choose(i,2)*length(which(temp==i))
    }
  }
  
  ## number of mutual nearest neighbors * 2
  mutual = length(which((nearest[nearest]-1:(2*n))==0))
  
  ## number of pairs that do not share any endpoint
  pairs = 2*n*(2*n-1)-2*share-2*2*n+mutual
  
  ## E(D12^2)  (share: C1; pairs: N(N-3)-2C1+2C2)
  p = n/(2*n-1)
  EDsq = n*p + share/2*p + pairs*(n-1)/4/(2*n-3)*p
  ED = n*p
  VD = EDsq - ED^2
  CovD = n*p + pairs*(n-1)/4/(2*n-3)*p - ED^2
  Z12 = (D12-ED)/sqrt(VD)
  Z21 = (D21-ED)/sqrt(VD)
  rho = CovD/VD
  x = min(Z12,Z21)
  # p1 = 2*pnorm(x)
  p2 =  1-pmvnorm(lower=rep(x,2),mean=rep(0,2),corr=matrix(c(1,rho,rho,1),2))[1]
  if (perm<=0){
    return(list(test.stat.D=min(D12,D21), test.stat.Z=x, pval.appr=p2))
  }else{
    stat = rep(0,perm)
    for (i in 1:perm){
      stat[i] = min(getD(distM, sample(1:N,n)))
    }
    p3 = length(which(stat<min(D12,D21)))/perm
    return(list(test.stat.D=min(D12,D21), test.stat.Z=x, pval.appr=p2, pval.perm=p3))
  }
}

# getD calculates D12 and D21

getD = function(distM, treated.index){
  N = dim(distM)[1]
  case.index = (1:N)[-treated.index]
  reorder = c(treated.index,case.index)
  distM2 = distM[reorder,reorder]
  diag(distM2) = max(distM)+10
  nearest = apply(distM2, 1, which.min)
  n = length(treated.index)
  D11 = length(which(nearest[1:n]<=n))
  D22 = length(which(nearest[(n+1):(2*n)]>n))
  D12 = n-D11 # number of nodes in sample 1 whose nearest neighbor is from sample 2
  D21 = n-D22
  return(c(D12,D21))
}

