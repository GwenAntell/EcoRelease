
# error handling: return NA values so won't cause dependent operations to fail

summarize_distrib <- function(distrib, denom=1, proportion=FALSE, include_singles=TRUE,
                              na.rm=FALSE, return_quarts=FALSE){
  if (include_singles==FALSE) { distrib <- distrib[distrib>1] }
  
  if (na.rm==TRUE){ distrib <- distrib[!is.na(distrib)] }
  
  if (length(distrib)==0){ 
    if (return_quarts==TRUE){ return(rep(NA, 5)) } else { return(rep(NA, 4)) } 
    } else {
    
    if (proportion==TRUE) {
      distrib <- distrib/denom
      if (any(distrib==1)){ distrib <- distrib-0.001 }
      distrib <- log10(distrib/(1-distrib)) # logit transformaiton
    }
    
    if (return_quarts==TRUE){
      output <- c(quantile(distrib, probs=c(0.25, 0.50, 0.75)), mean(distrib), sd(distrib))
      names(output) <- c('Q1','Q2','Q3','mean','sd')
    } else {
      output <- c(min(distrib), median(distrib), max(distrib), mean(distrib))
      names(output) <- c('min','med','max','mean')
    }
    return(output)
  }
}
