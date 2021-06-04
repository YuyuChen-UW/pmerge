##############################################################################################
##############################################################################################
##########Cauchy combination method#############################################################

# function inside the integral on LHS
f=function(t, epi, K){
  (K-1)*tan(pi*(1-epi+(K-1)*t-1/2))+tan(pi*(1-t-1/2))
}

#LHS-RHS
L_R=function(x ,epi, K){
  integrate(f, lower=x, upper=(epi/K),epi,K)$value-(epi/K-x)*f(x,epi,K)
}

# a_{F}

a_c=function(x, epi, K){
  H = (K-1)* tan(pi*((1-epi+(K-1)*x)-0.5)) + tan(pi*((1-x)-0.5))
  v = (1/pi)*atan(-(H)/K)+0.5
  return(v)
}


#' @title The Cauchy Merging Function
#' @description A function to merge arbitrarily dependent or independent p-values via the Cauchy combination method.
#' @param p numeric vector of p-values.
#' @param method method used to merge p-values, "A" or "I"; by default, method = "A".
#' @param epi The significance level of the hypothesis test; takes value in (0,0.5); by default epi = 0.1.
#' @param censor logical; if TRUE, p-values are left-censored at 0 and right-censored at 1; by default, censor = TRUE.
#' @details  The Cauchy combination method \insertCite{L}{pmerge} calculates the test statistic via
#' \deqn{ M_{ \mathcal C,K}(p_{1},\dots,p_{K}):=\mathcal C\left(\frac{1}{K}\sum_{k=1}^{K}\mathcal C^{-1}\left(p_{k}\right)\right).}
#' Here \eqn{\mathcal C} is the cumulative distribution function of the standard Cauchy random variable and \eqn{\mathcal C^{-1}} is the inverse of \eqn{\mathcal C}.
#' The hypothesis test is rejected if the test statistic is less than or equal to a threshold. The threshold is different for different methods [dependence="I" or "A"].
#' @details  "I" dependence [dependence="I"]: The threshold of the Cauchy combination method for independent p-values is (asymptotically) equivalent to the significance level; see \insertCite{L;textual}{pmerge}.
#' @details  "A" dependence [dependence="A"]: See \insertCite{C}{pmerge} for details about the threshold of the Cauchy combination method for arbitrarily dependent p-values.
#' @return The function returns the details of the hypothesis test including the method used, test result (i.e., reject or fail to reject), the significance level, the test statistic and the threshold.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{C}{pmerge}
#'
#' \insertRef{L}{pmerge}
#' @export

#The Cauchy combination method
pCauchy <- function(p, method="A",epi=0.1, censor = TRUE){
  K = length(p)

  if (censor == TRUE){
    p = p[ p>=0 & p<=1 ]
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  if (epi < 0 | epi > 0.5){
    warning("The significance level epi should be between 0 and 0.5")
  }

  val = tan(pi*(p-1/2))
  test.stat = atan(mean(val) )/pi+1/2

  if (method == "I"){
    threshold = epi

  } else if (method == "A"){
    #solve  x
    up = (epi/K)-1e-40
    low = 1e-40
    if (low>=up){
      stop("Bisection root searching: The lower bound is larger than the upper bound")
    }
    mid = (up+low)/2
    valmid = L_R(mid,epi = epi,K = K)
    while(abs(valmid)> 1e-4){
      if(valmid>0){
        up = mid
      }else{
        low = mid
      }
      mid = (up+low)/2
      valmid = L_R(mid,epi = epi,K = K)
      #print(valmid)
    }
    x = mid
    threshold = a_c(x = x,epi = epi,K = K)
  }

  #uniroot(f=L_R, lower = 1e-15, upper =  (epi/K)-1e-15, epi=0.1, K=K)

  significance.level = epi
  decision = ifelse(test.stat<=threshold,"Reject","Fail to reject")
  Results = c(method,decision,significance.level,test.stat,threshold)
  names(Results) = c("Method","Decision","Significance.level","Test.stat","Threshold")

  return(as.data.frame(Results))
}


#P<-runif(1000000,0,0.0001)
#pCauchy(P,method="A",epi=0.01)




