##############################################################################################
##############################################################################################
##########Hommel, Simes and grid harmonic#####################################################

#function for the grid harmonic method to solve ( use uniroot; gives the merged p-value)
gh <- function(epi,p){
  x=p/epi
  K=length(x)
  lk<-sum(1/c(1:K))
  val = sapply(c(1:K),function(i){ifelse(lk*x[i]<=1,1,0)})
  val = K*val/ceiling(K*lk*x)
  return(mean(val)-1)
}

#' @title The Simes Merging Function.
#' @description A function to merge arbitrarily dependent p-values via the Hommel method and the grid harmonic merging method,  or independent p-values via the Simes method.
#' @param p numeric vector of p-values.
#' @param method method to merge p-values, "H", "S" or "G"; by default, method = "H".
#' @param censor logical; if TRUE, p-values are left-censored at 0 and right-censored at 1; by default, censor = TRUE.
#' @param lower the lower bound for \code{\link{uniroot}} function; by default lower = 1e-10.
#' @details  Both "H" and "G" methods are valid for arbitrarily dependent p-values. The "S" method is valid for independent p-values.
#' @details "H" method [method="H"]: The Hommel method \insertCite{H}{pmerge} calculates the merged p-value by
#' \deqn{H_{K}(p_{1},\dots,p_{K})=\ell_{K}\wedge_{k=1}^{K}G_{k,K}(p_{1},\dots,p_{k}),}
#' where \eqn{p_{(k)}} is the \eqn{k}-th order statistic of the p-values, \eqn{\ell_{K}=\sum_{k=1}^{K}k^{-1}} and
#' \deqn{G_{k,K}(p_{1},\dots,p_{K})=\frac{K}{k}p_{(k)}\wedge 1.}
#' @details "S" method [method="S"]: The Simes method \insertCite{S}{pmerge} calculates the merged p-value by
#' \deqn{S_{K}(p_{1},\dots,p_{K})=\frac{1}{\ell_{K}}H_{K}(p_{1},\dots,p_{K}),}
#' where \eqn{H_{K}} is given in the details of "H" method.
#' @details "G" method [method="G"]: The grid harmonic method in \insertCite{VWW;textual}{pmerge}, which is an improvement of the Hommel method, produces a smaller merged p-value than the Hommel method; see \insertCite{VWW;textual}{pmerge} for more details.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @return The function returns the merged p-value.
#' @importFrom Rdpack reprompt
#' @note
#' Although the Simes method is developed for independent p-values, it is also valid for a wide range of dependency on p-values, e.g., MTP\eqn{_{2}} condition \insertCite{Sarkar}{pmerge}. However, it is not always valid for arbitrarily dependent p-values.
#' @references
#' \insertRef{H}{pmerge}
#'
#' \insertRef{Sarkar}{pmerge}
#'
#' \insertRef{S}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#' @export
pSimes <- function(p, method ="H", censor = TRUE, lower=10^(-10)){
  K = length(p)

  if (method == "G" & K > 1000000){
    warning("For large K: the calculation time for the grid harmonic method increases linearly as the number of p-values increases")
  }

  if (censor == TRUE){
    p[p<0]=0
    p[p>1]=1
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  if (method == "H"){
    pp = sum(1/c(1:K))*min((K/c(1:K))*sort(p),1)
  } else if (method == "S"){
    pp = min((K/c(1:K))*sort(p))
  } else if (method == "G"){
    # In the case of a p that is F=G <=>  Hk* = Hk and method = "G"
    # This case gives an error since gh(upper) < 0 hence same sign as gh(lower)
    # But in the paper of vovk et al. admissible it states:
    # "domination is strict if, in addition, F(p) < G(p) for at least one p"
       
    if (gh(epi = sum(1/c(1:K)) * min((K/c(1:K)) * sort(p), 1),p = p)<=0|gh(epi = lower,p = p)>0) {
      #message("Returning Hk instead of Hk*")
      pp = sum(1/c(1:K)) * min((K/c(1:K)) * sort(p), 1) # Hk
    }else{ # Hk* < Hk strictly smaller
      pp = uniroot(f = gh, lower = lower, upper = sum(1/c(1:K)) * min((K/c(1:K)) * sort(p), 1), p = p)$root #Hk*
    }
  }
  return(min(pp,1))
}

#test
#p=runif(10000)
#st = proc.time()
#phommel(p, method = "H")
#phommel(p, method = "S")
#phommel(p, method = "G")
#proc.time() - st
