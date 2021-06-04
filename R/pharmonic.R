##############################################################################################
##############################################################################################
##########  harmonic and harmonic* method under arbitrary dependence##########################

#' @title The Harmonic Mean Merging Function
#' @description A function to merge p-values via the harmonic mean and the harmonic* merging methods.
#' @param p numeric vector of p-values.
#' @param method method to merge p-values, "H1", "H2" or "H3"; by default, method = "H1".
#' @param censor logical; if TRUE, p-values are left-censored at 0 and right-censored at 1; by default, censor = TRUE.
#' @details  Both "H1" and "H2" methods are valid for dependent p-values and they require the number of p-values to be greater than 2. "H3" method is valid to merge independent p-values.
#' @details "H1" method [method="H1"]: The harmonic mean method for merging arbitrarily dependent p-values \insertCite{VW}{pmerge} calculates the merged p-value by
#' \deqn{F_{-1,K}(p_{1},\dots,p_{K})=b_{-1,K}M_{-1,K}(p_{1},\dots,p_{K}),}
#' where \eqn{M_{-1,K}(p_{1},\dots,p_{K})=\left(\frac{1}{K}\sum_{k=1}^{K}p_{k}^{-1}\right)^{-1}} is the harmonic mean of the p-values and \eqn{b_{-1,K}} is a constant; see \insertCite{VW;textual}{pmerge} for the calculation of \eqn{b_{-1,K}}.
#' @details "H2" method [method="H2"]: The harmonic\eqn{^{*}} merging method \insertCite{VWW}{pmerge}, which is an improvement of the harmonic mean method \insertCite{VW}{pmerge}, returns a smaller merged p-value than the harmonic mean method; see \insertCite{VWW;textual}{pmerge} for more details.
#' @details "H3" method [method="H3"]: See \insertCite{W;textual}{pmerge} about the harmonic mean method for merging independent p-values. The result is approximated by generalized central limit theorem; see \code{\link{pmean}} for details.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @importFrom Rdpack reprompt
#' @return The function returns the merged p-value.
#' @references
#' \insertRef{VW}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#'
#' \insertRef{W}{pmerge}
#' @export
pharmonic <- function(p, method = "H1", censor = TRUE){
  K = length(p)

  if (K == 2){
    stop("The number of p-values should be greater than 2")
  }

  if (censor == TRUE){
    p[p<0]=0
    p[p>1]=1
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  if (method == "H1"){
    pp = min(h(K)[1]*gmean(p,-1), 1)
  } else if (method == "H2"){
    n <- c(h(K)[2],rep(h(K)[3],K-1))
    p1=sort(p)^(-1)
    p2=n^(-1)
    m = sapply(c(1:K),function(i){sum(p2[1:i])/sum(p1[1:i])})
    #check
    #sp=sort(p)
    # m1=c()
    # for (i in 1:K){
    #   m1[i] = gmean(sp[1:i],-1)/gmean(n[1:i],-1)
    #  }
    # m-m1>=0.0001

    pp = min(m,ifelse(min(p) <=0, 0, 1))
  } else if (method == "H3"){
    pp = pmean(p=p, r=-1, dependence = "I")
  }
  return(pp)
}

#test
#p=runif(10^2)
#pharmonic(p, method = "H1")
#pharmonic(p, method = "H2")
