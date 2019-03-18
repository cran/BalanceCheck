
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

CrossMST = function(distM,treated.index,perm=0,k=1,discrete.correction=TRUE){
  n = length(treated.index)
  N = dim(distM)[1]
  case.index = (1:N)[-treated.index]
  reorder = c(treated.index,case.index)
  distM2 = distM[reorder,reorder]
  E = mstree(as.dist(distM2),k)
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
  m = N-n
  mu1 = k*n*(n-1)/N
  mu2 = k*m*(m-1)/N
  V1 = n*(n-1)*m*(m-1)/(N*(N-1)*(N-2)*(N-3))*((n-2)/(m-1)*(nEi+2*k*(N-1)-4*k^2*(N-1)^2/N)+k*(N-1)*(N-2*k)/N)
  V2 = n*(n-1)*m*(m-1)/(N*(N-1)*(N-2)*(N-3))*((m-2)/(n-1)*(nEi+2*k*(N-1)-4*k^2*(N-1)^2/N)+k*(N-1)*(N-2*k)/N)
  V12 = n*(n-1)*m*(m-1)/(N*(N-1)*(N-2)*(N-3))*(-nEi+k*(N-1)*(4*k*N-N-6*k)/N)

  temp = getR1R2(E,1:n)
  R1 = temp$R1
  R2 = temp$R2
  
  if (discrete.correction){
    R1 = R1-0.5
    R2 = R2-0.5
  }
  
  Z1 = (R1-mu1)/sqrt(V1)
  Z2 = (R2-mu2)/sqrt(V2)
  rho = V12/sqrt(V1*V2)
  x = max(Z1,Z2)
  p2 =  1-pmvnorm(upper=rep(x,2),mean=rep(0,2),corr=matrix(c(1,rho,rho,1),2))[1]
  if (perm<=0){
    return(list(test.stat.Z=x, pval.appr=p2))
  }else{
    Rmat = Zmat = matrix(0,perm,2)
    for (i in 1:perm){
      temp2 = getR1R2(E,sample(1:N,n))
      Rmat[i,] = c(temp2$R1,temp2$R2)
      Zmat[i,] = c(Rmat[i,]-c(mu1,mu2))/c(sqrt(V1),sqrt(V2))
    }
    stat = apply(Zmat,1,max)
    p3 = length(which(stat>=x))/perm
    return(list(test.stat.Z=x, pval.appr=p2,  pval.perm=p3))
  }
}



## function CrossNN calculates the test statistic and p-value for the new test
## if perm>0, also calculate permutation p-value

CrossNN = function(distM,treated.index,perm=0,k=1,discrete.correction=TRUE){
  n = length(treated.index)
  N = dim(distM)[1]
  
  if (k==1){
    case.index = (1:N)[-treated.index]
    reorder = c(treated.index,case.index)
    distM2 = distM[reorder,reorder]
    diag(distM2) = max(distM)+10
    nearest = apply(distM2, 1, which.min)
    n = length(treated.index)
    D11 = length(which(nearest[1:n]<=n))
    D22 = length(which(nearest[(n+1):N]>n))
    m = N-n
    
    ### pick one arrow, pick another arrow, 2n*(2n-1) possibilities
    
    ## the number of possibilities that the two arrows share the endpoint
    temp = as.vector(table(nearest))
    a = max(temp)
    share = 0
    if (a>1){
      for (i in 2:a){
        share = share + choose(i,2)*length(which(temp==i))
      }
    }
    
    ## number of mutual nearest neighbors * 2
    mutual = length(which((nearest[nearest]-1:N)==0))
  }else{
    temp = getDk(distM, treated.index, k)
    A = temp$A
    D11 = temp$D11
    D22 = temp$D22
    temp2 = table(A)
    id = as.numeric(row.names(temp2))
    deg = rep(0,N)
    deg[id] = temp2
    share = (sum(deg^2)-sum(deg))/2
    count = 0
    for (i in 1:N){
      ids = A[i,]
      count = count + length(which(A[ids,]==i))
    }
    mutual = count
  }

  ## C1: mutual/2; C2: share
  
  ED11 = k*n*(n-1)/(N-1)
  ED22 = k*m*(m-1)/(N-1)
  VD11 = n*m*(n-1)*(m-1)/(N*(N-1)*(N-2)*(N-3))*(k*N+mutual+(n-2)/(m-1)*(share*2+k*N-k^2*N)-2*k^2*N/(N-1))
  VD22 = n*m*(n-1)*(m-1)/(N*(N-1)*(N-2)*(N-3))*(k*N+mutual+(m-2)/(n-1)*(share*2+k*N-k^2*N) -2*k^2*N/(N-1))
  CovD = n*m*(n-1)*(m-1)/(N*(N-1)*(N-2)*(N-3))*(mutual-2*share+k^2*N*(N-3)/(N-1))
  
  if (discrete.correction){
    D11 = D11-0.5
    D22 = D22-0.5
  }
  
  Z1 = (D11-ED11)/sqrt(VD11)
  Z2 = (D22-ED22)/sqrt(VD22)
  rho = CovD/sqrt(VD11*VD22)
  x = max(Z1,Z2)
  # p1 = 2*pnorm(x)
  p2 =  1-pmvnorm(upper=rep(x,2),mean=rep(0,2),corr=matrix(c(1,rho,rho,1),2))[1]
  if (perm<=0){
    return(list(test.stat.Z=x, pval.appr=p2))
  }else{
    Dmat = Zmat = matrix(0,perm,2)
    for (i in 1:perm){
      if (k==1){
        temp2 = getD(distM, sample(1:N,n))
        Dmat[i,] = temp2
      }else{
        temp2 = getDk(distM, sample(1:N,n), k)
        Dmat[i,] = c(temp2$D11, temp2$D22)
      }
      Zmat[i,] = (Dmat[i,]-c(ED11,ED22))/c(sqrt(VD11),sqrt(VD22))
    }
    stat = apply(Zmat,1,max)
    p3 = length(which(stat>=x))/perm
    return(list(test.stat.Z=x, pval.appr=p2, pval.perm=p3))
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
  D22 = length(which(nearest[(n+1):N]>n))
  return(c(D11,D22))
}

getDk = function(distM, treated.index, k){
  N = dim(distM)[1]
  case.index = (1:N)[-treated.index]
  reorder = c(treated.index,case.index)
  distM2 = distM[reorder,reorder]
  diag(distM2) = max(distM)+10
  
  A = matrix(0,N,k)
  for (i in 1:N){
    A[i,] = (sort(distM2[i,1:N], index.return=T)$ix)[1:k]
  }
  
  n = length(treated.index)
  D11 = length(which(A[1:n,]<=n))
  D22 = length(which(A[(n+1):N,]>n))
  return(list(A=A,D11=D11,D22=D22))
}

