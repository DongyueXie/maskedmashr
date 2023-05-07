
#'@title Run mash
#'@param adjust how to adjust the covariance matrix, either prior or lb or const. lb is for lower bound
#'@importFrom mashr get_significant_results
#'@importFrom mashr cov_pca
#'@importFrom mashr mash
#'@importFrom mashr mash_1by1
#'@export
mash_wrapper = function(Bhat,P=NULL,npc=5,Shat=NULL,
                        nu = ncol(Bhat)+1,verbose=FALSE,
                        adjust = "lb",
                        U.true = NULL,
                        U.ed=NULL,
                        n.ed=NULL,
                        adj.const = 0.1){



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

    if(is.null(U.ed)){

      m.1by1 = mash_1by1(data)
      strong = get_significant_results(m.1by1)
      n.ed = length(strong)
      if(length(strong)==0){
        strong=NULL
        n.ed=N
      }
      U.pca = cov_pca(data,npc,strong)
      U.ed = mashr:::bovy_wrapper(data,U.pca,subset=strong)
    }
    U.est = U.ed$Ulist
    ##############
    ##############
    # U.est  = lapply(U.est,function(x){
    #   x - diag(1/sqrt(N),R)
    # })
    #############
    #############
    if(!is.null(adjust)){

      if(adjust == 'prior'){
        U.est = lapply(1:length(U.est),function(k){
          nk = n.ed*U.ed$pi[k]
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
          nk = n.ed*U.ed$pi[k]
          if(nk>1){
            U.est[[k]] + 2*diag(sqrt(2/nk),R)
            #U.est[[k]] + 2*diag(sqrt(2/nk*1/(diag(solve(U.est[[k]]+diag(R))))^2))
            #temp = eigen(U.est[[k]])
            #temp$vectors%*%(diag(pmax(2/sqrt(nk),temp$values)))%*%t(temp$vectors)
          }else{
            NULL
          }
        })
      }

      if(adjust == "const"){
        U.est = lapply(1:length(U.est),function(k){
          U.est[[k]] + diag(adj.const,R)
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
