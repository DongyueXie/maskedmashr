#'@title Function for running masked mash on a single condition.
#'@param Z matrix of Z scores
#'@param P matrix of p values, can be NULL, and will be calculated from Z scores
#'@param r the condition that we are interested in for FDR control
#'@param p.thresh threshold of p-values for masking. It take values in [0,0.5]. p<=thresh and p>=(1-thresh) will be masked.
#'@param npc number of data driven covariance matrices for mash
#'@param adjust method for adjusting diagonal of estimated covariance matrices, "lb" or "prior".
#'@param strong_for_md Whether use strong effects for estimating data driven cov matrices
#'@param strong.thresh for deciding strong effects.
#'@importFrom mashr mash_set_data
#'@importFrom mashr cov_canonical
#'@export

maskedmash_condition_wrapper = function(Z,P=NULL,
                              r,
                              p.thresh=0.5,
                              npc=5,
                              adjust = 'lb',
                              adj.const = 1e-3,
                              verbose=FALSE,
                              strong_for_md = TRUE,
                              strong.thresh = 0.2,
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
  P.temp = pmax(P.temp,1e-12)
  P.temp = pmin(P.temp,1-1e-12)

  Z = sign(Z)*qnorm(1-P.temp/2)
  if(is.null(p.thresh)){
    p.thresh = 0.5
  }
  thresh = qnorm(1-p.thresh/2)
  Z = mash_set_data(Z)
  U.c = cov_canonical(Z)
  if(strong_for_md){
    if(verbose){
      cat("Finding strong effects for deconvolution using masked data...")
      cat('\n')
    }
    ## get strong signals for estimating data-driven matrices
    Z.mask = t(apply(Z$Bhat,1,function(x){
      x.alt = mask.func(x)
      is.mask = is.mask.z(x,thresh)
      is.mask[-r] = FALSE
      is.mask*sign(x)*pmax(abs(x),abs(x.alt)) + (1-is.mask)*x
    }))
    strong = maskedmash_get_strong(Z.mask,thresh,strong.thresh = strong.thresh)
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
  Z.comb = mask.Z(Z$Bhat,thresh,r)
  U.est = masked.md(Z$Bhat,Z.comb,strong=strong,thresh=thresh,
                    usepointmass = usepointmass,U.canon=U.c,
                    U.data=NULL,npc=npc,adjust=adjust,verbose=verbose,adj.const=adj.const)$U.est.adj
  out = masked.mash(Z$Bhat,Z.comb,thresh=thresh,U.canon = U.c,U.data = U.est,verbose=verbose,return_post_weights=return_post_weights)
  out$p.thresh = p.thresh
  out$P = P
  out$r = r
  t1 = Sys.time()
  out$run_time = difftime(t1,t0)
  out
}
