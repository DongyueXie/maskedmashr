#'@title simulate Z scores
#'@param N total number of samples
#'@param Ulist a list of correlation matrices
#'@param Pi prior weights
#'@param prior normal, t, uniform
#'@param mean.range range of alternative E(z), from -mean.range to mean.range
#'@importFrom mvtnorm rmvnorm
#'@importFrom mvtnorm rmvt
#'@importFrom mvtnorm qmvt
#'@export
simDataI.ult = function(N,Ulist,Pi=NULL,prior="t",t.df=10,mean.range = 4){

  K = length(Ulist)
  R = ncol(Ulist[[1]])
  if(is.null(Pi)){Pi = rep(1/K,K)}
  nsamp = round(N*Pi)
  if(sum(nsamp)<N){
    nsamp[K] = nsamp[K] + N - sum(nsamp)
  }
  B = c()
  Bhat = c()

  if(prior=='normal'){

    for(k in 1:K){
      Bk = rmvnorm(nsamp[k],rep(0,R),(mean.range/qnorm(0.99))^2*Ulist[[k]])
      B = rbind(B,Bk)
    }

  }else if(prior == 't'){

    for(k in 1:K){
      Bk = rmvt(nsamp[k],sigma=(mean.range/qmvt(0.99,df=t.df,tail="lower.tail")$quantile)^2*Ulist[[k]],df=t.df)
      B = rbind(B,Bk)
    }

  }else if(prior == 'uniform'){

    for(k in 1:K){
      Bk = rmvu(nsamp[k],Ulist[[k]])
      B = rbind(B,Bk)
    }
    non0.idx = which(B!=0)
    B[non0.idx] = B[non0.idx]*2*mean.range-mean.range
  }else if(prior == 'half.uniform'){
    for(k in 1:K){
      Bk = rmvu(nsamp[k],Ulist[[k]])
      B = rbind(B,Bk)
    }
    B = B*mean.range
  }else{
    stop('prior should be one of: normal, t, uniform, half.uniform')
  }


  N = sum(nsamp)
  Bhat = matrix(rnorm(N*R,B,sd=1),ncol = R)
  P = 2*(1-pnorm(abs(Bhat)))
  Shat = matrix(1,nrow=N,ncol=R)

  list(B=B,P=P,Bhat=Bhat,Shat=Shat,input = list(N=N,Ulist=Ulist,Pi=Pi,prior=prior,t.df=t.df,mean.range = mean.range))


}



#'@title generating multivariate uniform(0,1) distribution with correlation U
rmvu = function(n,U){
  non0.idx = which(diag(U)!=0)
  if(length(non0.idx)==ncol(U)){
    U = cov2cor(U)
    X = rmvnorm(n,rep(0,ncol(U)),U)
    X = apply(X,2,pnorm)
  }else if(length(non0.idx)==0){
    X = matrix(0,nrow=n,ncol=ncol(U))
  }else{
    U.non0 = U[non0.idx,non0.idx,drop=F]
    U.non0 = cov2cor(U.non0)
    X = matrix(0,nrow=n,ncol=ncol(U))
    x = rmvnorm(n,rep(0,ncol(U.non0)),U.non0)
    x = apply(x,2,pnorm)
    X[,non0.idx] = x
  }
  X
}
