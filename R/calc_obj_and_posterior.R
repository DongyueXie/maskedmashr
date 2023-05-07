#'@title calculate posterior weights
#'@param Z.comb augmented z score matrices
#'@param lik.list output of `calc_relative_lik_matrix` from mashr
#'@param Pi prior weights
#'@return posterior weights
#'

calc_post_weights = function(Z.comb,lik.list,Pi){

  N = length(Z.comb)
  post_weights = lapply(1:N,function(i){
    # lik.mat is m*K
    lik.mat = exp(lik.list[[i]]$loglik_matrix+lik.list[[i]]$lfactors)
    # gamma_i is m*K
    gamma_i = t(t(lik.mat*Z.comb[[i]]$z.det.jacob)*Pi)
    gamma_i/sum(gamma_i)
  })

  post_weights

}

#'@title evaluate objective function
calc_obj_maskedmd = function(lik.list,Z.comb,Pi,N,R,multiplier = NULL){
  if(is.null(multiplier)){
    multiplier = rep(1,length(Pi))
  }
  l = sum(log(unlist(lapply(1:N,function(i){
    lm_i = lik.list[[i]]
    det_i = Z.comb[[i]]$z.det.jacob
    sum(det_i*exp(mashr:::compute_vloglik_from_matrix_and_pi(Pi,lm_i, matrix(1,nrow=length(det_i),ncol=R))))
  }))))
  idx = which(Pi!=0)
  l + sum((multiplier[idx]-1)*log(Pi[idx]))
}


#'@title calculate posterior summaries, given data and prior parameters.
#'@importFrom ashr compute_lfsr
calc_post_summary = function(Z.comb,Ulist,Pi,post_weights,pi_thresh){

  N = length(Z.comb)
  R = nrow(Ulist[[1]])

  I_R = diag(R)
  which.comp = which(Pi > pi_thresh)
  Ulist = Ulist[which.comp]
  Pi = Pi[which.comp]
  post_weights = lapply(post_weights,function(x){x[,which.comp,drop=FALSE]})

  K = length(Ulist)

  result = vector(mode='list',length = 6)
  names(result) = c('PosteriorMean','PosteriorSD','lfdr','lfsr','NegativeProb','PosteriorWeights')

  #post_weights = vector(mode='list',length = N)
  post_mean = matrix(nrow=N,ncol=R)
  post_cov = vector(mode='list',length = N)
  lfdr = matrix(nrow=N,ncol=R)
  negativeProb = matrix(nrow = N,ncol = R)

  # pre-calculation of U(U+I)^{-1}
  post_cov_byk = lapply(Ulist,function(x){x%*%solve(x+I_R)})
  #post_cov_byk = lapply(Ulist,function(x){solve(ginv(x)+I_R)})
  # K*R indicator matrix : indicate 0 diagonal value
  post_cov_0diag = lapply(post_cov_byk,function(x){diag(x)==0})
  post_cov_0diag = matrix(unlist(post_cov_0diag),ncol=R,byrow = TRUE)


  for(i in 1:N){
    #For each z score vector, generate all possible combinations
    # Z.comb_i = enumerate.z(Z[i,],thresh)
    Z.comb_i  = Z.comb[[i]]
    Ji = nrow(Z.comb_i$z.comb)

    # 1. calculate posterior weights
    # gamma_i = matrix(nrow=Ji,ncol = K)
    # for(j in 1:Ji){
    #   for(k in 1:K){
    #     gamma_i[j,k] = pi[k]*dmvnorm(Z.comb_i$z.comb[j,],sigma = I_R+Ulist[[k]])*Z.comb_i$z.det.jacob[j]
    #   }
    # }
    # gamma_i = gamma_i/sum(gamma_i)
    # post_weights[[i]] = gamma_i
    gamma_i = post_weights[[i]]

    # 2. calculate posterior mean and variance
    mu.tilde_i = 0
    U.tilde_i = 0
    negativeProb_i = matrix(0,nrow=K,ncol = R)
    #positiveProb_i = matrix(0,nrow=K,ncol = R)
    for(k in 1:K){
      for(j in 1:Ji){
        mu.tilde_ijk = post_cov_byk[[k]]%*%Z.comb_i$z.comb[j,]
        mu.tilde_i = mu.tilde_i + gamma_i[j,k]*mu.tilde_ijk
        U.tilde_i = U.tilde_i + gamma_i[j,k]*(post_cov_byk[[k]]+tcrossprod(mu.tilde_ijk))
        negativeProb_i[k,] = negativeProb_i[k,] + gamma_i[j,k]*pnorm(-mu.tilde_ijk/sqrt(diag(post_cov_byk[[k]])))
        #positiveProb_i[k,] = positiveProb_i[k,] + gamma_i[j,k]*(1-pnorm(-mu.tilde_ijk/sqrt(diag(post_cov_byk[[k]]))))
      }
    }
    post_mean[i,] = mu.tilde_i
    post_cov[[i]] = U.tilde_i - tcrossprod(mu.tilde_i)

    # 3. calculate lfdr

    lfdr[i,] = colSums(colSums(gamma_i)*post_cov_0diag)

    # 4. calculate negativeProb
    negativeProb[i,] = colSums(negativeProb_i*(1-post_cov_0diag),na.rm = TRUE)
  }

  lfsr = compute_lfsr(negativeProb,lfdr)

  result$PosteriorMean = post_mean
  result$PosteriorSD = sqrt(matrix(unlist(lapply(post_cov,diag)),ncol=R,byrow = TRUE))
  result$lfdr = lfdr
  result$lfsr = lfsr
  result$NegativeProb = negativeProb



  result

}
