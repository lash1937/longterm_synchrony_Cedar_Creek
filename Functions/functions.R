

#' Smith and Wilson's Evenness Index (Evar)
#' 
#' This function computes the Evar index of evenness from a vector of species abundances
#' for a given plot or community unit
#'
#' @param abundances 
#'
#' @return Numeric index between 0 and 1
#' @export
#'
#' @examples
#' 
#' abundances <- rpois(10, lambda = 10)
#' Evar(abundances)
#' 
Evar <- function(abundances){
  
  if(sum(abundances == 0) > 0){
    warning(
      "Removing species with zero abundance"
    )
  }
  abund_sub <- abundances[abundances != 0]
  
  if(length(abund_sub) <= 1){
    warning(
      "Returning NA: At least two species must be present to compute evenness"
    )
    return(NA)
  }
  # log the abundances
  l_n <- log(abund_sub)
  n <- length(l_n)
  
  
  # compute variance
  vln <- var(l_n) * (n - 1) / n
  
  # transform to constrain to [0, 1]
  return(
    1 - (2 / pi) * atan(vln)
  )
  
}


#' Compute Box-Cox Transformation for a given lambda
#'
#' @param x Vector of strictly positive data.
#' @param lambda Chosen lambda based on profile likelihood plot
#'
#' @return Box-Cox transformed x.
#' @export
#'
#' @examples
#' 
boxcox_transform <- function(x, lambda){
  
  if(sum(x <= 0) > 0){
    stop("Box-Cox Transform only appropriate for strictly positive data.
         Some observations <= 0.")
  }
  if(lambda == 0){
    return(log(x))
  } else{
    return((x^lambda - 1) / lambda)
  }
  
}

aictable<-function(X,m){
  rnames<-row.names(X)
  AICc<-X$AIC+2*X$df*(X$df+1)/(m-X$df-1)     #small-sample correction
  logL<-X$df-X$AIC/2                         #Log-likelihood
  tab<-data.frame(X[,1],logL,AICc)           #Remove AIC column; add logL and AICc
  colnames(tab)[1]<-c("Params")              #Rename "df" column
  row.names(tab)<-rnames
  tab<-tab[order(tab$AICc),]                 #Sort by ascending AICc value
  deltaAICc<-tab$AICc-min(tab$AICc)          #Delta AICc
  weight<-exp(-deltaAICc/2)/sum(exp(-deltaAICc/2))  #Weights
  cumwt<-weight                              #Column for cumulative weight
  for(i in 2:dim(X)[1]){
    cumwt[i]<-cumwt[i-1]+cumwt[i]              #Accumulate weight from the top
  }
  tab<-data.frame(tab,deltaAICc,weight,cumwt)
  tab<-round(tab,4)
  tab
}
