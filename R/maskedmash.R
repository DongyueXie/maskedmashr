#'@title Wrapper function for running masked mash
#'@param Z matrix of Z scores
#'@param P matrix of p values, can be NULL, and will be calculated from Z scores
#'@param p.thresh threshold of p-values for masking. It take values in [0,0.5]. p<=thresh and p>=(1-thresh) will be masked.
#'@param npc number of data driven covariance matrices for mash
#'@param adjust method for adjusting diagonal of estimated covariance matrices, "lb" or "prior".
#'@param strong_for_md Whether use strong effects for estimating data driven cov matrices
#'@param qval.thresh for deciding strong effects.
#'@importFrom mashr mash_set_data
#'@importFrom mashr cov_canonical
#'@export

maskedmash_wrapper = function(Z,P=NULL,
                              p.thresh=0.5,
                              npc=5,
                              adjust = 'lb',
                              verbose=FALSE,
                              strong_for_md = TRUE,
                              qval.thresh = 0.3,
                              usepointmass = TRUE,
                              return_post_weights = FALSE){
  t0 = Sys.time()

  if(is.null(Z)){
    stop('Z must be provided')
  }
  if(is.null(P)){
    P = 2*(1-pnorm(abs(Z)))
  }

  R = ncol(Z)
  N = nrow(Z)

  P.temp = P
  P.temp = pmax(P.temp,1e-10)
  P.temp = pmin(P.temp,1-1e-10)

  Z = sign(Z)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)
  data = mash_set_data(Z)
  U.c = cov_canonical(data)
  #m.1by1 = mash_1by1(data)
  #strong = get_significant_results(m.1by1,0.2)
  if(strong_for_md){
    if(verbose){
      cat("Finding strong effects for deconvolution using AdaPT with masked data...")
      cat('\n')
    }
    ## get strong signals for estimating data-driven matrices
    strong = maskedmash_get_strong(P,p.thresh,qval_thresh = qval.thresh)
    if(verbose){
      cat(paste("Found",length(strong),"strong effects",sep=" "))
      cat('\n')
    }
    if(length(strong)==0){
      strong=NULL
    }
  }else{
    strong=NULL
  }

  #U.pca = cov_pca(data,npc,strong)
  #browser()
  U.est = masked.md(data,strong=strong,thresh=thresh,
                    usepointmass = usepointmass,U.canon=U.c,
                    U.data=NULL,npc=npc,adjust=adjust,verbose=verbose)$U.est.adj
  out = masked.mash(data,thresh=thresh,U.canon = U.c,U.data = U.est,verbose=verbose,return_post_weights=return_post_weights)
  out$p.thresh = p.thresh
  out$P = P
  t1 = Sys.time()
  out$run_time = difftime(t1,t0)
  out
}

#'@title get strong effects index
#'@importFrom matrixStats rowMins
#'@importFrom adaptMT adapt_glm
#'@export
maskedmash_get_strong = function(P,p.thresh,qval_thresh=0.2){
  qv = matrix(nrow=nrow(P),ncol=ncol(P))
  for(r in 1:ncol(P)){
    qv[,r] = adapt_glm(x=data.frame(x=rep(1,nrow(P))),
                    pvals = P[,r],
                    pi_formulas = "x",
                    mu_formulas = "x",
                    nfits=1,
                    s0=rep(p.thresh,nrow(P)),
                    alphas = 0,
                    verbose = list(print = FALSE, fit = FALSE, ms = FALSE))$qvals
  }
  return(which(rowMins(qv)<=qval_thresh))
}

#'@title masked mash
#'@param data object from mash_set_data
#'@param thresh NULL or a number, |z-score| larger than thresh will be masked.
#'@param Pi prior weights; either fixed or as init value
#'@param U.canon a list of canonical prior cov matrices; fixed.
#'@param U.data a list of data-driven prior cov marrices; either fixed or as initialization.
#'@param fixg whether fix the prior weights and covariance matrices.
#'@param usepointmass whether to include a point mass at 0
#'@param pi_thresh threshold below which mixture components are ignored
#'@param normalizeU whether normalize U such that the largest diag element is 1
#'@param algorithm.version Rcpp or R for evaluating likelihood
#'@param optmethod 'mixSQP' or 'EM'
#'@param verbose TRUE to print progress
#'@param return_post_weights whether return the full posterior weights. Could be large
#'@param control a list of control parameters for SQP
#'@param prior nullbiased or uniform
#'@return a list of Posterior mean, sd, lfsr, lfdr, negativeProb.
#'@export
masked.mash = function(data,
                       thresh=NULL,
                       Pi=NULL,
                       U.canon=NULL,
                       U.data,
                       fixg=FALSE,
                       usepointmass=TRUE,
                       #U.update = 'none',
                       pi_thresh = 1e-5,
                       gridmult = sqrt(2),
                       grid = NULL,
                       normalizeU = TRUE,
                       algorithm.version = 'Rcpp',
                       optmethod = 'mixSQP',
                       prior = 'nullbiased',
                       max_iter=1000,
                       tol=1e-5,
                       verbose=TRUE,
                       printevery = 100,
                       return_post_weights = FALSE,
                       control=list()){
  Z = data$Bhat
  N = nrow(Z)
  R = ncol(Z)
  I_R = diag(R)

  if(is.null(thresh)){
    # mask all a-scores
    thresh = qnorm(0.75)
  }
  Z.comb = mask.Z(Z,thresh)

  # if prior is fixed, then directly calculate posterior summaries

  ## posterior weights: of dimension N*J_i*K. Use list, each list element is a J_i*K matrix
  ##
  if(fixg){

    Ulist = c(U.canon,U.data)
    if(normalizeU){Ulist = mashr:::normalize_Ulist(Ulist)}
    if(!is.null(grid)){
      xUlist = mashr:::expand_cov(Ulist,grid,usepointmass)
    }else{
      if(usepointmass){
        xUlist = c(list(null=matrix(0,nrow = R,ncol = R)),Ulist)
      }else{
        xUlist = Ulist
      }
    }

    if(is.null(Pi)){Pi = rep(1/length(xUlist),length(xUlist))}
    lik.list = lapply(Z.comb,function(x){
      temp = mashr:::calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = algorithm.version)
      temp
    })

    # calculate posterior weights
    post_weights = calc_post_weights(Z.comb,lik.list,Pi)
    #browser()
    # calculate posteriors
    #result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
    #result$PriorCov = Ulist
    #result$PosteriorWeights = post_weights
    loglik = calc_obj(lik.list,Z.comb,Pi,N,R)

  }else{

    ############## estimate weights ####################
    if(is.null(grid)){
      grid = mashr:::autoselect_grid(data,gridmult)
    }
    if(verbose){
      cat("Estimating prior weights...")
      cat("\n")
    }
    out = estimate_pi(Z.comb,U.canon,U.data,grid,usepointmass,
                      Pi,max_iter,tol,pi_thresh,normalizeU,
                      algorithm.version,optmethod,verbose,printevery,control,prior)
    Pi = out$pi
    xUlist = out$xUlist
    post_weights = out$post_weights
    loglik = out$loglik

  }


  # get posterior summaries
  if(verbose){
    cat("Calculating posterior distributions...")
    cat("\n")
  }
  result = calc_post_summary(Z.comb,xUlist,Pi,post_weights,pi_thresh)


  effect_names = rownames(data$Bhat)
  condition_names = colnames(data$Bhat)
  for (i in names(result)) {
    if (length(dim(result[[i]])) == 2) {
      colnames(result[[i]]) = condition_names
      rownames(result[[i]]) = effect_names
    }
  }

  out = list()
  out$result = result
  out$maskedProp = sum(log(unlist(lapply(Z.comb,function(x){nrow(x$z.comb)})),2))/(N*R)
  out$loglik = loglik
  out$fitted_g = list(pi=Pi,Ulist = c(U.canon,U.data),grid = grid, usepointmass = usepointmass)
  if(return_post_weights){
    out$posterior_weights = post_weights
  }
  out

}



#'@title Estimate prior weights using masked data, given covariance matrices.
estimate_pi = function(Z.comb,U.canon,U.data,grid,usepointmass,Pi,max_iter,tol,
                       pi_thresh,normalizeU,algorithm.version,optmethod,verbose,printevery,control,prior){

  N = length(Z.comb)
  R = ncol(U.data[[1]])
  n_grid = length(grid)

  Ulist = c(U.canon,U.data)
  if(normalizeU){Ulist = mashr:::normalize_Ulist(Ulist)}
  K = length(Ulist)
  xUlist = mashr:::expand_cov(Ulist,grid,usepointmass)
  P = length(xUlist) # the order of P is : null, (l=1,k=1:K), (l=2,k=1:K),...
  if(is.null(Pi)){Pi = rep(1/P,P)}
  # marginal likelihood
  ## need to include L
  ### marginal likelihood list: a list of length N, each element is a matrix of dim Ji*P.
  ### posterior weights: a list of length N, each element is a matrix of dim Ji*P
  ### prior weights pi is a length P vector

  # calculate likelihood
  lik.list = lapply(Z.comb,function(x){
    temp = mashr:::calc_relative_lik_matrix(mash_set_data(x$z.comb),xUlist,algorithm.version = algorithm.version)
    temp
  })

  if(optmethod=='mixSQP'){

    # Get the likelihood matrix, of dimension (sum_i Ji) * K
    lik_mat = lapply(1:N,function(i){
      colSums(Z.comb[[i]]$z.det.jacob * exp(lik.list[[i]]$loglik_matrix+lik.list[[i]]$lfactors))
    })
    lik_mat = do.call(rbind,lik_mat)
    prior = mashr:::set_prior(ncol(lik_mat),prior)
    Pi = mashr:::optimize_pi(lik_mat,prior=prior,optmethod=optmethod,control=control)

  }

  if(optmethod=='EM'){

    if(prior == 'uniform'){
      multiplier = rep(1,P)
    }
    if(prior == 'nullbiased'){
      multiplier = c(10,rep(1,P-1))
    }

    loglik.obj = -Inf

    for(iter in 1:max_iter){

      # evaluate objective function
      loglik.obj[iter+1] = calc_obj_maskedmd(lik.list,Z.comb,Pi,N,R,multiplier)

      ## check convergence
      if((loglik.obj[iter+1]-loglik.obj[iter])/N<tol){
        break
      }

      # update gamma_ijkl
      post_weights = calc_post_weights(Z.comb,lik.list,Pi)


      ## update pi
      Pi = lapply(post_weights,colSums)
      nk = Reduce('+',Pi) + multiplier - 1
      Pi = nk/sum(nk)

      if(verbose){
        if(iter%%printevery==0){
          cat(sprintf("Done %1.0f iterations, loglikelihood at %.2f",iter,loglik.obj[iter+1]))
          cat("\n")
        }
      }

    }
  }



  #result = calc_post_summary(Z.comb,xUlist,pi,post_weights,pi_thresh)
  result = list()
  result$post_weights = calc_post_weights(Z.comb,lik.list,Pi)
  result$xUlist = xUlist
  result$pi = Pi
  result$loglik = calc_obj_maskedmd(lik.list,Z.comb,Pi,N,R)

  result

}

