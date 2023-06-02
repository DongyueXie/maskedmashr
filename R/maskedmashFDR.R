#'@title masked mash for FDR control
#'@description order the tests from most to least significant based on lfsr from masked.mash, and remove least significant ones until FDP<alpha
#'@param obj fitted maskedmash object, or a list of P, result$lfsr, p.thresh
#'@param alpha target FDR level
#'@return a list of rej.set indx and fdp_hat
#'@export
maskedmashFDR = function(obj,alpha = 0.05){

  P = obj$P
  N = nrow(P)
  R = ncol(P)
  score = obj$result$lfsr
  rej.region = which(is.mask.p(P,obj$p.thresh))
  P = P[rej.region]
  score = score[rej.region]

  # rank the tests from most significant to least significant
  # and also as initial rejection set
  rej.set = rej.region[order(score,decreasing = FALSE)]


  fdp.t = fdp.hat(rej.set,obj$P)
  fdp.tr = fdp.t
  rej.idx = length(rej.set)
  rej.idx.lb = 1
  rej.idx.rb = length(rej.set)


  if(fdp.t>alpha){

    while(fdp.t>alpha | (fdp.t<=alpha&fdp.tr<=alpha)){

      if(rej.idx.lb==rej.idx.rb){
        rej.idx = NULL
        break
      }

      rej.idx = floor((rej.idx.lb+rej.idx.rb)/2)
      fdp.t = fdp.hat(rej.set[1:rej.idx],obj$P)
      fdp.tr = fdp.hat(rej.set[1:(rej.idx+1)],obj$P)

      if(fdp.t>alpha){
        rej.idx.rb = rej.idx
      }

      if((fdp.t<=alpha&fdp.tr<=alpha)){
        rej.idx.lb = rej.idx
      }

    }

  }

  if(is.null(rej.idx)){
    list(rej.set=NULL,fdp.hat=0,alpha=alpha)
  }else{
    rej.set=rej.set[1:rej.idx]
    rej.set=rej.set[obj$P[rej.set]<0.5]
    list(rej.set=rej.set,fdp.hat=fdp.t,alpha=alpha)
  }

}

maskedmashFDR_condition = function(obj,alpha = 0.05){
  P = obj$P[,obj$r]
  #N = nrow(P)
  #R = ncol(P)
  score = obj$result$lfsr[,obj$r]
  rej.region = which(is.mask.p(P,obj$p.thresh))
  P = P[rej.region]
  score = score[rej.region]

  # rank the tests from most significant to least significant
  # and also as initial rejection set
  rej.set = rej.region[order(score,decreasing = FALSE)]


  fdp.t = fdp.hat(rej.set,obj$P[,obj$r])
  fdp.tr = fdp.t
  rej.idx = length(rej.set)
  rej.idx.lb = 1
  rej.idx.rb = length(rej.set)


  if(fdp.t>alpha){

    while(fdp.t>alpha | (fdp.t<=alpha&fdp.tr<=alpha)){

      if(rej.idx.lb==rej.idx.rb){
        rej.idx = NULL
        break
      }

      rej.idx = floor((rej.idx.lb+rej.idx.rb)/2)
      fdp.t = fdp.hat(rej.set[1:rej.idx],obj$P[,obj$r])
      fdp.tr = fdp.hat(rej.set[1:(rej.idx+1)],obj$P[,obj$r])

      if(fdp.t>alpha){
        rej.idx.rb = rej.idx
      }

      if((fdp.t<=alpha&fdp.tr<=alpha)){
        rej.idx.lb = rej.idx
      }

    }

  }

  if(is.null(rej.idx)){
    list(rej.set=NULL,fdp.hat=0,alpha=alpha)
  }else{
    rej.set=rej.set[1:rej.idx]
    rej.set=rej.set[(obj$P[,obj$r])[rej.set]<0.5]
    list(rej.set=rej.set,fdp.hat=fdp.t,alpha=alpha)
  }

}

#'@title Estimate of FDP as in adapt paper, without threshold version
fdp.hat = function(rej.set,P){

  R = sum(P[rej.set]<=1/2)
  A = sum(P[rej.set]>1/2)
  (1+A)/(max(R,1))

}

