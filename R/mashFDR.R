#'@title mash for FDR control
#'@description order the tests from most to least significant based on lfsr, and reject least significant ones until FDP<alpha
#'@param obj fitted mash object
#'@param alpha target FDR level

mashFDR = function(obj,alpha = 0.05){

  lfsr = obj$result$lfsr
  # rank the tests from most significant to least significant
  # and also as initial rejection set
  if(mean(lfsr)<=alpha){
    rej.set = 1:length(lfsr)
  }else{
    l = order(lfsr,decreasing = FALSE)
    lfsr.ordered  =lfsr[l]
    lfsr.cmean = cumsum(lfsr.ordered)/(1:length(lfsr.ordered))
    rej.set = l[1:(which(lfsr.cmean>alpha)[1]-1)]
  }

  return(rej.set)

}
