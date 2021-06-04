#' pmerge
#'
#' @description The pmerge-package contains functions to merge p-values that are independent or arbitrarily dependent.
#' @details
#' \enumerate{
#' \item
#' If p-values are independent, methods for merging independent p-values should be used.
#' \item
#' If the dependence structure of p-values is unknown, methods for merging arbitrarily dependent p-values should be used.
#' Arbitrarily dependence of p-values means that p-values can have any dependence structure such as independence, perfectly positive dependence and so on.
#' Regardless of the dependence structure of p-values, methods for merging arbitrarily dependent p-values always produce a valid p-value in the sense that the probability of making a Type-I error is below the significance level; see, e.g., \insertCite{VW;textual}{pmerge} for details.
#' \item
#' This package contains the following functions/methods to merge p-values:
#' \itemize{
#' \item The Generalized Mean Merging Function (\code{\link{pmean}}): It contains the generalized mean methods for independent p-values \insertCite{C}{pmerge} and arbitrarily dependent p-values (\insertCite{VW;textual}{pmerge} and \insertCite{VWW;textual}{pmerge}).
#' \item The Order Statistics Merging Function (\code{\link{porder}}): It contains the order statistics merging method \insertCite{VWW}{pmerge} for arbitrarily dependent p-values.
#' \item The Harmonic Mean Merging Function (\code{\link{pharmonic}}): It contains the harmonic mean method \insertCite{W}{pmerge} for independent p-values, the harmonic mean method \insertCite{VW}{pmerge} and the harmonic* merging method \insertCite{VWW}{pmerge} for arbitrarily dependent p-values.
#' \item The Simes Merging Function (\code{\link{pSimes}}): It contains the Simes method \insertCite{S}{pmerge} for independent p-values, the Hommel method \insertCite{H}{pmerge} and the grid harmonic merging method \insertCite{VWW}{pmerge} for arbitrarily dependent p-values.
#' \item The Cauchy Merging Function (\code{\link{pCauchy}}):  It contains the Cauchy combination methods for independent p-values \insertCite{L}{pmerge} and arbitrarily dependent p-values \insertCite{C}{pmerge}.
#' }
#' }
#' @references
#' \insertRef{C}{pmerge}
#'
#' \insertRef{H}{pmerge}
#'
#' \insertRef{L}{pmerge}
#'
#' \insertRef{S}{pmerge}
#'
#' \insertRef{VW}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#'
#' \insertRef{W}{pmerge}
#' @docType package
#' @name pmerge-package
NULL
