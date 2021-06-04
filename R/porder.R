##############################################################################################
##############################################################################################
##########Order statistics method#############################################################
#' @title The Order Statistics Merging Function
#' @description A function to merge arbitrarily dependent p-values via order statistics merging method.
#' @param p numeric vector of p-values.
#' @param k the \eqn{k}-th order statistic, takes value \eqn{1,2,\dots,K} where \eqn{K} is the number of p-values.
#' @param censor logical; if TRUE, p-values are left-censored at 0 and right-censored at 1; by default, censor = TRUE.
#' @details The order statistics merging method combines p-values using order statistics and it includes the famous Bonferroni method [k=1, method="O1"].
#' @details It is first proposed by \insertCite{R;textual}{pmerge} to merge arbitrarily dependent p-values by
#' \deqn{G_{k,K}(p_{1},\dots,p_{K})=\frac{K}{k}p_{(k)}\wedge 1,}
#' where \eqn{p_{(k)}} is the \eqn{k}-th order statistic of the p-values.
#' @details This function uses an improved version \insertCite{VWW}{pmerge} of the above order statistics merging method. The improved method returns a smaller merged p-value than the order statistics merging method \insertCite{R}{pmerge}. The merged p-value is calculated by
#' \deqn{G^{*}_{k,K}(p_{1},\dots,p_{K})=G_{k,K}(p_{1},\dots,p_{K})\wedge 1_{\{p_{(1)}>0\}}.}
#' where \eqn{p_{(1)}} is the first order statistic of the p-values.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @return The function returns the merged p-value.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{R}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#' @export
porder <- function(p, k, censor = TRUE){
  K = length(p)

  if (censor == TRUE){
    p[p<0]=0
    p[p>1]=1
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  pp = min(min((K/k)*sort(p)[k],1),ifelse(min(p)>0,1,0))

  return(pp)
}

#test
#p=runif(10000)
#porder(p,k=5, method = "O1")
#porder(p,k=5, method = "O2")
