#'@title masked deconvolution: estimate prior covariance matrices using multivariate deconvolution and masked data
#'@description This function estimates prior weights and variances using masked z scores. If strong
#'is provided, then U.canon should be set to NULL and usepointmass set to FALSE;
#'@param data object from `mash_set_data` function
#'@param strong index of samples that have "strong" effects.
#'@param thresh NULL or a number, |z-score| larger than thresh or |z|<=s(t) will be masked. when thresh<=Phi^{-1}(0.75), all masked.
#'@param Pi prior weights; either fixed or use as init value
#'@param U.canon a list of canonical prior cov matrices; fixed.
#'@param U.data a list of data-driven prior cov matrices; If null, will use pca on masked z score for initialization.
#'@param algorithm.version either Rcpp or R for calculating likelihood
#'@param adjust adjust diagonal of estimated U
#'@param nu prior df of U. If NULL, nu = R+1, for adjust = 'prior'
#'@importFrom mashr cov_pca
masked.md = function(data,
                     strong=NULL,
                     thresh=NULL,
                     Pi=NULL,
                     U.canon=NULL,
                     U.data = NULL,
                     npc = 5,
                     usepointmass=FALSE,
                     algorithm.version = 'Rcpp',
                     adjust = 'lb',
                     nu = NULL,
                     pi_thresh = 1e-8,
                     max_iter=1000,
                     tol=1e-5,
                     verbose=TRUE,
                     printevery = 50){

  Z = data$Bhat
  if(is.null(strong)){
    strong = 1:nrow(Z)
  }else{
    usepointmass = FALSE
    U.canon = NULL
  }
  Z = Z[strong,]
  N = nrow(Z)
  R = ncol(Z)
  I_R = diag(R)

  if(is.null(thresh)){
    # mask all z-scores
    thresh = qnorm(0.75)
  }

  Z.comb = mask.Z(Z,thresh)
  if(is.null(U.data)){
    U.data = cov_pca(mash_set_data(do.call(rbind, lapply(Z.comb, function(x) x$z.mask))),npc)
  }
  Ulist = c(U.canon,U.data)
  if(usepointmass){
    Ulist = c(list(null=matrix(0,nrow = R,ncol = R)),Ulist)
  }
  K = length(Ulist)

  if(is.null(Pi)){Pi = rep(1/K,K)}
  loglik.obj = -Inf

  if(verbose){
    cat("Running deconvolution")
    cat("\n")
  }
  for(iter in 1:max_iter){
    # calculate marginal likelihood
    # a list of length N, each a J*K log likelihood matrix
    lik.list = lapply(Z.comb,function(x){
      temp = mashr:::calc_relative_lik_matrix(mash_set_data(x$z.comb),Ulist,algorithm.version = algorithm.version)
      temp
    })

    # evaluate objective function
    # compute_vloglik_from_matrix_and_pi returns a m*1 loglik vector(after summing over k)
    loglik.obj[iter+1] = calc_obj_maskedmd(lik.list,Z.comb,Pi,N,R)
    ## check convergence
    if(abs(loglik.obj[iter+1]-loglik.obj[iter])/N<tol){
      break
    }

    # update gamma_ijk
    post_weights = calc_post_weights(Z.comb,lik.list,Pi)

    ## update pi
    Pi = lapply(post_weights,colSums)
    Pi = Reduce('+',Pi)/N

    # update U.data
    ## eigenvectors of sum_ij(gamma_ijk S_ij)
    ## let gamma.tilde_ij = sum_l gamma_ijkl w_l/(1+w_l)

    for(k in (usepointmass+length(U.canon)+1):K){

      S.tilde.k = lapply(1:N,function(i){
        crossprod(sqrt(post_weights[[i]][,k])*Z.comb[[i]]$z.comb)
      })
      #browser()
      S.tilde.k = Reduce("+",S.tilde.k)/(N*Pi[k])

      temp = eigen(S.tilde.k)
      Ulist[[k]] = temp$vectors%*%diag(pmax(temp$values-1,0))%*%t(temp$vectors)
    }

    if(verbose){
      if(iter%%printevery==0){
        cat(sprintf("Done %1.0f iterations, loglikelihood at %.2f",iter,loglik.obj[iter+1]))
        cat("\n")
      }
    }



  }

  ## now we have estimated prior weights and covariance
  ## calculate posterior summaries

  #result = calc_post_summary(Z.comb,xUlist,pi,post_weights)

  U.est = Ulist[(usepointmass+length(U.canon)+1):K]
  pi.est = Pi[(usepointmass+length(U.canon)+1):K]

  if(!is.null(adjust)){
    if(adjust == 'prior'){
      if(is.null(nu)){
        nu = R+1
      }
      U.est.adj = lapply(1:length(U.est),function(k){
        nk = N*pi.est[k]
        if(nk>R){
          s2.hat = R*nu/(nk+nu)/sum(diag(solve(nk*(U.est[[k]]+I_R))))
          (nk/(nk+nu-R-1))*U.est[[k]] + ((R+1-nu+s2.hat)/(nk+nu-R-1))*I_R
        }else{
          NULL
        }
      })

    }

    if(adjust == "lb"){
      U.est.adj = lapply(1:length(U.est),function(k){
        nk = N*pi.est[k]
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
  }

  U.est.adj = U.est.adj[lengths(U.est.adj) != 0]
  result = list()
  result$U.est = U.est
  result$U.est.adj = U.est.adj
  result$pi = Pi
  result$loglik = loglik.obj
  result$post_weights = post_weights

  result
}
