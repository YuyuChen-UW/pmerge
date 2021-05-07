#install.packages("stabledist")
library(stabledist)
##############################################################################################
##############################################################################################
########## Generalized mean function##########################################################

# The generalized mean
# A function to calculate the generalized mean of a vector.
# p numeric vector of values.
# r the exponent of the generalized mean function.
# return the generalized mean of the vector.


gmean<-function(p,r){
  if (r == 0) {
    m <- exp(mean(log(p)))
  } else if(r == -Inf){
    m <- min(p)
  } else if(r == Inf){
    m <- max(p)
  } else {
    m <- mean(p^(r))^(1/r)
  }
  return(m)
}


##############################################################################################
##############################################################################################
########## Multipliers for the averaging method under arbitrary dependence####################

#Inputs:
# K: the number of p-values
# r: the exponent of the averaging function

#r=0: multiplier for geometric mean
#asymptotic result is used for K>=20
g<-function(K){

  if (K < 20){
  y = uniroot(function(y){
   K-(K^(2))/(y+K)-log(y+1)
  },lower=0.000000001,upper=10^300)$root
  } else {
    y = 0
  }

  c = 1/(y+K)
  a <- c()
  if (K < 20){
    a = (1/c)*exp(-(K-1)*(1-K*c))
  } else {
    a = exp(1)
  }

  return(a)
}
#test: ratio g(k)/exp(k) (Table 2 of "Combining p-values via averaging")
#for (i in 2:25){
#  x=g(i)
#  print(i)
#  print(x/exp(1))}

#r=-1 & K>2 : multiplier for harmonic mean (first element of the output)
#asymptotic result is not used due to the slow convergence rate

h<-function(K){
  y = uniroot(function(y){
    y^(2)-K*((y+1)*log(y+1)-y)
  },lower=0.00001,upper=10^100)$root

  c=1/(y+K)
  d = 1-(K-1)*c
  a=((y+K)^(2))/((y+1)*K)
  return(c(a,c,d))
}


#test: ratio h(k)/log(k) (Table 3 of "Combining p-values via averaging")
#for (i in 3:100000){
#  x=h(i)[1]
#  print(i)
#  print(x/log(i))}
#h(20000000)[1]/log(20000000)

#multiplier for r<1/(K-1) & r!= 0 or -1 & K>2 (r is restricted to r<-1)
#equation to solve in the function n(r,K) ( replace c=1/(y+K) )
LR<-function(y,r,K){
  val = (K-1)*((y+K)/(y+1))^(-r)+(y+K)^(-r)-K*(((y+K)/(y+1))^(-r-1)-(y+K)^(-r-1))/((r+1)*(y/(y+K)))
}
# for K \geq 3000000, asymptotic result is used
n<-function(r,K){
  if (K <=3000000){
    if (r<0){

      if (r > (-5) & r <=(-2)){
        up = 10^(K/(10000*r*r)) + K
      } else if (r<=-5){
        up = 10^(K/(1000*r*r)) + K
      } else if (r>= (-2) & r <= (-1.5)) {
        up = 10^(K/(10000*r*r)) + K^(2)
      } else if (r>= (-1.5) & r < 0) {
        up = 10^(K/(10000*r*r)) + K^(10)
      }

      low = 1e-7

      mid = (up+low)/2
      valmid = LR(mid,r,K)
      eps = up-low

      while(abs(valmid) > 0.001 &  eps > 10^(-5)){
        #&  eps > 10^(-11)
        if(valmid > 0){
          up = mid
        }else{
          low = mid
        }
        mid = (up+low)/2
        valmid = LR(mid,r,K)
        eps = up-low
      }
      y = mid
    }else {
      y=uniroot(LR,lower=10,upper=10^11,K=K,r=r)[1]
      y=as.numeric(y)
    }
    c=1/(y+K)
    d = 1-(K-1)*c
    a = 1/gmean(c(c,rep(d,times=(K-1))),r)
  } else if (K >3000000){
    a =  (r/(r+1))*K^(1+(1/r))
  }


  return(a)
}

#test
#for (i in 3:10000){
#  x=n(-10,i)
#  print(i)
#  print(x)}
#n(-5,20000000)

#test: Check the asymptotic results for r<-1
#K=1000000
#R<-seq(-10,-1.1,by=0.1) #r values
#result<-c()
#for (i in 1: length(R) ){
#result[i]<-n(R[i],K)/(  (R[i]/(R[i]+1))*K^(1+(1/R[i]))  )
#}
#result[abs(result-1) > 0.0001]
#result

#the multiplier a_{r,K} for arbitrary dependence
amean<-function(r,K){
  if (K == 2 & r< 1) {
    a <- 2
  } else if (r == -Inf){
    a <- K
  } else if (r == Inf){
    a <- 1
  } else if (r >= (1/(K-1)) & r < Inf ){
    a <- min(r+1,K)^(1/r)
  } else if (r==0){
    a <- g(K)
  } else if (r==-1){
    a <- h(K)[1]
  } else {
    a <-n(r,K)
  }
  return(a)
}


#par(mfrow=c(2,3))

#R<-seq(-10,10,by=0.1)
#K=10000
#A=c()
#for (i in 1:length(R)){
#  A[i]=amean(R[i],K)
#  print(i)
#}

#plot(R,A,main="K=10000",ylab=expression(a[r]),xlab="r")

#amean(-0.2,100000)

##############################################################################################
##############################################################################################
########## The averaging method under arbitrary dependence and independence###################




#(1): r is rounded to the closest integer for arbitrary dependence assumption, large K and -2 < r < 0
#     (avoid multiplier convergence issues and computational issues for r close to -1 or 0)
#(2): if r < 0, zero p-values are omitted.
#(3): asymptotic result is used for r < -1 and large K



#' @title The Generalized Mean Merging Function
#' @description A function to merge dependent p-values or independent p-values via the generalized mean method.
#' @param p numeric vector of p-values.
#' @param r the exponent of the generalized mean function in the set \eqn{[-10,\infty]\cup\{-\infty\}};
#' @param dependence the dependence assumption of p-values, "A" or "I".
#' @param subset logical; if TRUE, p-values less than 0 or greater than 1 are removed; by default, subset=FALSE.
#' @details  The generalized mean merging method combines p-values through the generalized mean function
#' \deqn{M_{r,K}(p_{1},\dots,p_{K})=\left(\frac{1}{K}\sum_{k=1}^{K}p_{k}^{r}\right)^{\frac{1}{r}},}
#' including the limiting cases:
#' \deqn{M_{-\infty,K}(p_{1},\dots,p_{K})=\min\{p_{1},\dots,p_{K}\};}
#' \deqn{M_{0,K}(p_{1},\dots,p_{K})=\left(\prod_{k=1}^{K}p_{k}\right)^{\frac{1}{K}};}
#' \deqn{M_{\infty,K}(p_{1},\dots,p_{K})=\max\{p_{1},\dots,p_{K}\}.}
#' @details Some special cases of the generalized mean merging method are:
#'
#'  (1) The Bonferroni method [\eqn{r=-\infty}, method="A"];
#'
#'
#'  (2) The harmonic averaging method [\eqn{r=-1}, method="A"] \insertCite{VW}{pmerge};
#'
#'  (3) The Fisher method [\eqn{r=0} and method="I"].

#' @details  "A" dependence [dependence="A"]: The method is valid for dependent p-values. The merged p-value is calculated by
#' \deqn{F_{r,K}(p_{1},\dots,p_{K})=b_{r,K}M_{r,K}(p_{1},\dots,p_{K}),}
#' where \eqn{M_{r,K}} is the generalized mean of the p-values with exponent \eqn{r}, and \eqn{b_{r,K}} is a constant; the constant \eqn{b_{r,K}} can be found in \insertCite{VW;textual}{pmerge} and \insertCite{VWW;textual}{pmerge}.
#' @details  "I" dependence [dependence="I"]: The method is valid for independent p-values. The merged p-value is approximated by the classic central limit theorem or the generalized central limit theorem; see \insertCite{C;textual}{pmerge} for more details. For exponent less than \eqn{0}, large number of p-values is needed to ensure the accuracy of approximation.
#' @details  When the number of p-values is large, if \eqn{r\in[-2,1]} and method="A", \eqn{r} is rounded to the closest integer to prevent computational issues and convergence issues; see \insertCite{VW;textual}{pmerge} and \insertCite{VWW;textual}{pmerge}.
#' @details  Zero p-values are omitted if r is less than zero
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @import stabledist
#' @import stats
#' @importFrom Rdpack reprompt
#' @return The function returns the merged p-value.
#' @references
#' \insertRef{C}{pmerge}
#'
#' \insertRef{VW}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#' @export
pmean<-function( p, r = 0, dependence = "A", subset = FALSE){
  K=length(p)

  if (subset == TRUE){
    p = p[ p>=0 & p<=1 ]
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  #  if (r != -Inf & r-round(r, digits=1)!=0 & r != Inf){
  #    warning("r is rounded to the first digit")
  #    r=round(r,digits=1)
  #  }


  if (r > (-2) & r< (1) & dependence == "A" & K >= 10000){
    warning("r is rounded to the closest integer if r is between -2 and 1 and dependence = 'A' for large number of p-values")
    r=round(r)
  }

  if (r < -10 & r > (- Inf)){
    stop("r should be greater than -10 or equal to -Inf")
  }

  if (r < 0 & min(p) <= 0){
    warning("Zero p-values are omitted if r is less than or equal to zero")
    p = p[p>0]
  }

  if (dependence == "A") {
    if (K == 2 & r< 1) {
      pp = gmean(p,r) * 2
    } else {
      pp = gmean(p,r) * amean(r,K)
    }
  } else if (dependence == "I") {
    if (r == 0){
      pp = 1 - pchisq(q = -2*K*log(gmean(p,r)),df=2*K)
    } else if (r == - Inf){
      pp = 1 - (1 - gmean(p,r))^(K)
    } else if (r == Inf){
      pp = gmean(p,r)^(K)
    } else if (r > 0 & r < Inf){
      mu = (r+1)^(-1)
      sig = (r^(2))*((1+2*r)^(-1))*((1+r)^(-2))
      pp = pnorm((gmean(p,r)^(r)-mu)*(sqrt(sig/K)))
    } else if (r > -Inf & r<0){

      alpha = (-1/r)
      if (r > (-1/2)){
        ca = (K*((alpha/(alpha-2))-(alpha/(alpha-1))^(2)))^(1/2)
        bk = (K*alpha)/(alpha-1)
      }
      else if (r == (-1/2)){
        ca = (K*log(K))^(0.5)
        bk = (K*alpha)/(alpha-1)
      }
      else if (r < (-1/2) & r>-1){
        ca = (K^(1/alpha))*(gamma(1-alpha)*cos(pi*alpha/2))^(1/alpha)
        bk = (K*alpha)/(alpha-1)
      } else if (r == -1) {
        ca = K*pi/2
        bk = ((K^(2))*pi/2)*integrate(f = function(x){ sin(2*x/(K*pi))*alpha*x^(-alpha-1)}, lower = 1, upper = Inf)$value
      } else {
        ca = (K^(1/alpha))*(gamma(1-alpha)*cos(pi*alpha/2))^(1/alpha)
        bk = 0
      }

      if (r >=(-1/2)){
        pp = 1- pnorm ((K*gmean(p,r)^(r)-bk)/ca)
      } else {
        pp = 1- pstable ((K*gmean(p,r)^(r)-bk)/ca,  alpha , beta = 1 , gamma = 1, delta = 0, pm = 1)
      }

    }
  }

  return(min(pp,1))
}



#test

#p=runif(10^4)
#p=c(p,0,-1)

#meanp(p, r=-2, dependence = "A",subset=TRUE)
#meanp(p, r=-2, dependence = "I",subset=TRUE)


##############################################################################################
##############################################################################################
##########  harmonic and harmonic* method under arbitrary dependence##########################

#(1): if any of the p-values does not satisfy 0 <= p <=1, a warning is given
#(2): number of p-values needs to be greater than 2 (for harmonic*)



#' @title The Harmonic Mean Merging Function
#' @description A function to merge dependent p-values via the harmonic and the harmonic* merging methods.
#' @param p numeric vector of p-values.
#' @param method method to merge p-values, "H1" or "H2"; by default, method="H1".
#' @param subset logical; if TRUE, p-values less than 0 or greater than 1 are removed; by default, subset=FALSE.
#' @details  Both "H1" and "H2" methods are valid for dependent p-values and they require the number of p-values to be greater than 2.
#' @details "H1" method [method="H1"]: The harmonic mean method \insertCite{VW}{pmerge} calculates the merged p-value by
#' \deqn{F_{-1,K}(p_{1},\dots,p_{K})=b_{-1,K}M_{-1,K}(p_{1},\dots,p_{K}),}
#' where \eqn{M_{-1,K}(p_{1},\dots,p_{K})=\left(\frac{1}{K}\sum_{k=1}^{K}p_{k}^{-1}\right)^{-1}} is the harmonic mean of the p-values and \eqn{b_{-1,K}} is a constant; see \insertCite{VW;textual}{pmerge} for the calculation of \eqn{b_{-1,K}}.
#' @details "H2" method [method="H2"]: The harmonic\eqn{^{*}} merging method \insertCite{VWW}{pmerge}, which is a modification of the harmonic averaging method \insertCite{VW}{pmerge}, returns a smaller merged p-value than the harmonic averaging method; see \insertCite{VWW;textual}{pmerge} for more details.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @importFrom Rdpack reprompt
#' @return The function returns the merged p-value.
#' @references
#' \insertRef{VW}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#' @export
pharmonic <- function(p, method = "H1", subset = FALSE){

  K = length(p)

  if (K == 2){
    stop("The number of p-values should be greater than 2")
  }

  if (subset == TRUE){
    p = p[ p>=0 & p<=1 ]
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }


  m <- c()
  n <- c(h(K)[2],rep(h(K)[3],K-1))

  for (i in 1:K){
    m[i] = gmean(sort(p)[1:i],-1)/gmean(n[1:i],-1)
  }

  if (method == "H1"){
    pp = min( h(K)[1]*gmean(p,-1)  , 1 )
  } else if (method == "H2"){
    pp = min(m,ifelse(min(p) <=0, 0, 1))
  }
  return(pp)
}

#test
#p=runif(10^2)
#harmonicp(p, method = "H1")
#harmonicp(p, method = "H2")

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

#' @title The Hommel Merging Function.
#' @description A function to merge dependent p-values via the Hommel method and the grid harmonic merging method,  or independent p-values via the Simes method.
#' @param p numeric vector of p-values.
#' @param method method to merge p-values, "H", "S" or "G"; by default, method="H".
#' @param subset logical; if TRUE, p-values less than 0 or greater than 1 are removed; by default, subset=FALSE.
#' @details  Both "H" and "G" methods are valid for dependent p-values. The "S" method is valid for independent p-values.
#' @details "H" method [method="H"]: The Hommel method \insertCite{H}{pmerge} calculates the merged p-value by
#' \deqn{H_{K}(p_{1},\dots,p_{K})=l_{K}\wedge_{k=1}^{K}G_{k,K}(p_{1},\dots,p_{k}),}
#' where \eqn{p_{(k)}} is the \eqn{k}-th order statistic of the p-values, \eqn{l_{K}=\sum_{k=1}^{K}k^{-1}} and
#' \deqn{G_{k,K}(p_{1},\dots,p_{K})=\frac{K}{k}p_{(k)}\wedge 1.}
#' @details "S" method [method="S"]: The Simes method \insertCite{S}{pmerge} calculates the merged p-value by
#' \deqn{S_{K}(p_{1},\dots,p_{K})=\frac{1}{l_{K}}H_{K}(p_{1},\dots,p_{K}),}
#' where \eqn{H_{K}} is given in the details of "H" method.
#' @details "G" method [method="G"]: The grid harmonic method in \insertCite{VWW;textual}{pmerge}, which is a modification of the Hommel method, produces a smaller merged p-value than the Hommel method; see \insertCite{VWW;textual}{pmerge} for more details.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @return The function returns the merged p-value.
#' @importFrom Rdpack reprompt
#' @note
#' Although the Simes method is developed for independent p-values, it is also valid for a wide range of dependency on p-values, e.g., MTP\eqn{_{2}} condition \insertCite{Sarkar}{pmerge}. However, it is not always valid for dependent p-values.
#' @references
#' \insertRef{H}{pmerge}
#'
#' \insertRef{Sarkar}{pmerge}
#'
#' \insertRef{S}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#' @export
phommel <- function(p, method ="H", subset = FALSE){
  K = length(p)

  if (method == "G" & K > 1000000){
    warning("For large K: the calculation time for the grid harmonic method increases linearly as the number of p-values increases")
  }

  if (subset == TRUE){
    p = p[ p>=0 & p<=1 ]
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
    pp=uniroot(f=gh,lower=0.00001,upper=sum(1/c(1:K))*min((K/c(1:K))*sort(p),1),p=p)$root
  }

  return(min(pp,1))

}

#test
#p=runif(10000)
#st = proc.time()
#hommelp(p, method = "H")
#hommelp(p, method = "S")
#hommelp(p, method = "G")
#proc.time() - st

##############################################################################################
##############################################################################################
##########Order statistics method#############################################################

#' @title The Order Statistics Merging Function
#' @description A function to merge dependent p-values via order statistics merging method.
#' @param p numeric vector of p-values.
#' @param k the \eqn{k}-th order statistic, takes value \eqn{1,2,\dots,K} where \eqn{K} is the number of p-values.
#' @param method method used to merge p-values, "O1" or "O2"; by default, method="O1".
#' @param subset logical; if TRUE, p-values less than 0 or greater than 1 are removed; by default, subset=FALSE.
#' @details The order statistics merging method combines p-values using order statistics and it includes the famous Bonferroni method [k=1, method="O1"].
#' @details  Both "O1" and "O2" methods are valid for dependent p-values.
#' @details "O1" method [method="O1"]: The order statistics class of merging methods proposed by \insertCite{R;textual}{pmerge} calculates the merged p-value by
#' \deqn{G_{k,K}(p_{1},\dots,p_{K})=\frac{K}{k}p_{(k)}\wedge 1,}
#' where \eqn{p_{(k)}} is the \eqn{k}-th order statistic of the p-values.
#' @details "O2" method [method="O2"]: It is a modification \insertCite{VWW}{pmerge} of the order statistics merging method \insertCite{R}{pmerge} by a zero-one adjustment. It returns a smaller merged p-value than the order statistics merging method \insertCite{R}{pmerge}. The merged p-value is calculated by
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
porder <- function(p, k, method="O1", subset = FALSE){
  K = length(p)

  if (subset == TRUE){
    p = p[ p>=0 & p<=1 ]
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  if  (method == "O1"){
    pp = min((K/k)*sort(p)[k],1)
  } else if (method== "O2"){
    pp = min(min((K/k)*sort(p)[k],1),ifelse(min(p)>0,1,0))
  }
  return(pp)
}

#test
#p=runif(10000)
#orderp(p,k=5, method = "O1")
#orderp(p,k=5, method = "O2")
