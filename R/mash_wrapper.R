
#'@title Run mash
#'@param adjust how to adjust the covariance matrix, either prior or lb or const. lb is for lower bound
#'@importFrom mashr get_significant_results
#'@importFrom mashr cov_pca
#'@importFrom mashr mash
#'@importFrom mashr mash_1by1
#'@export
mash_wrapper = function(Bhat,P=NULL,
                        npc=5,Shat=NULL,
                        nu = ncol(Bhat)+1,
                        verbose=FALSE,
                        adjust = "lb",
                        md.method = 'teem',
                        strong.thresh = 0.2,
                        U.true = NULL,
                        U.data=NULL,
                        n.data=NULL,
                        adj.const = 1e-3){



  if(!is.null(P)){
    Z = sign(Bhat)*qnorm(1-P/2)
  }else{
    Z = Bhat
  }

  R = ncol(Z)
  N = nrow(Z)

  data = mash_set_data(Z,Shat)
  U.c = cov_canonical(data)

  if(is.null(U.true)){

    if(is.null(U.data)){

      m.1by1 = mash_1by1(data)
      strong = get_significant_results(m.1by1,thresh = strong.thresh)
      n.data = length(strong)
      if(length(strong)==0){
        strong=NULL
        n.data=N
      }
      U.pca = cov_pca(data,npc,strong)
      if(md.method=='ed'){
        U.data = mashr:::bovy_wrapper(data,U.pca,subset=strong)
        U.est = U.data$Ulist
        pi.est = pi.est
      }else if(md.method=='teem'){
        U.data = mashr:::teem_wrapper(data,U.pca,subset=strong)
        U.est = U.data$U
        U.est = lapply(seq(dim(U.est)[3]), function(x) U.est[ , , x])
        pi.est = U.data$w
      }

    }

    ##############
    # U.est  = lapply(U.est,function(x){
    #   x - diag(1/sqrt(N),R)
    # })
    #############
    #############
    if(!is.null(adjust)){

      if(adjust == 'prior'){
        U.est = lapply(1:length(U.est),function(k){
          nk = n.data*pi.est[k]
          if(nk>R){
            s2.hat = R*nu/(nk+nu)/sum(diag(solve(nk*(U.est[[k]]+diag(R)))))
            (nk/(nk+nu-R-1))*U.est[[k]] + ((R+1-nu+s2.hat)/(nk+nu-R-1))*diag(R)
          }else{
            NULL
          }

        })
      }

      if(adjust == "lb"){
        U.est = lapply(1:length(U.est),function(k){
          nk = n.data*pi.est[k]
          if(nk>1){
            U.est[[k]] + 2*diag(sqrt(2/nk),R)
            #U.est[[k]] + 2*diag(sqrt(2/nk*1/(diag(solve(U.est[[k]]+diag(R))))^2))
            #temp = eigen(U.est[[k]])
            #temp$vectors%*%(diag(pmax(2/sqrt(nk),temp$values)))%*%t(temp$vectors)
          }else{
            U.est[[k]] + 2*diag(sqrt(2/n.data),R)
          }
        })
      }

      if(adjust == "const"){
        U.est = lapply(1:length(U.est),function(k){
          cov2cor(U.est[[k]]) + diag(adj.const,R)
        })
      }

    }

    U.data = U.est[lengths(U.est) != 0]

  }else{
    U.data=U.true
  }

  out = mash(data,c(U.c,U.data),verbose=verbose)
  out

}
