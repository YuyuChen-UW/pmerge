##############################################################################################
##############################################################################################
########## Generalized mean function##########################################################

# The generalized mean
# A function to calculate the generalized mean of a vector.
# p: numeric vector of values.
# r: the exponent of the generalized mean function.
# return the generalized mean of the vector.

gmean <- function(p,r){
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
# K: the number of p-values
# r: the exponent of the averaging function

# multiplier for geometric mean
# asymptotic result is used for K>=20
g <- function(K){
  if (K < 20){
    y = uniroot(function(y){
      K-(K^(2))/(y+K)-log(y+1)
    }, lower = 1e-9, upper=1e+300)$root
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
h <- function(K){
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
#h(2000000000)[1]/log(2000000000)

#multiplier for r<1/(K-1) & r!= 0 or -1 & K>2
#equation to solve in the function n(r,K) ( replace c = 1/(y+K) )
LR <- function(y,r,K){
  val = (K-1)*((y+K)/(y+1))^(-r)+(y+K)^(-r)-K*(((y+K)/(y+1))^(-r-1)-(y+K)^(-r-1))/((r+1)*(y/(y+K)))
}
# exact multiplier
n0 <- function(r,K){
  if (r<0){

    if (r > (-5) & r <=(-2)){
      up = 10^(K/(10000*r*r)) + K
    } else if (r<=-5){
      up = 10^(K/(1000*r*r)) + K
    } else if (r>= (-2) & r <= (-1.5)) {
      up = 10^(K/(10000*r*r)) + K^(2)
    } else if (r>= (-1.5) & r < 0) {
      up = 10^(K/(10000000*r*r)) + K^(10)
      #up = 10^(K/(10000*r*r)) + K^(10)
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
  } else {
    y=uniroot(LR,lower=10,upper=10^11,K=K,r=r)[1]
    y=as.numeric(y)
  }

  c=1/(y+K)
  d = 1-(K-1)*c
  a = 1/gmean(c(c,rep(d,times=(K-1))),r)
  return(a)
}
# asymptotic result is used
n <- function(r,K){
  if (r > -1.5 & r < -1 & K >= 6000000){
    a = (r/(r+1))*K^(1+(1/r))
  } else if (r > -1 & r < 1 & K >= 6000000){
    a = (r+1)^(1/r)
  } else if (r <= -1.5 & K>= 1000000 ){
    a = (r/(r+1))*K^(1+(1/r))
  } else if (r >=1 & K>= 1000000){
    a = (r+1)^(1/r)
  } else {
    a = n0(r,K)
  }
  return(a)
}
#test
#for (i in 3:10000){
#  x=n(-10,i)
#  print(i)
#  print(x)}
#n(-0.1,100000000)

#the multiplier a_{r,K} for arbitrary dependence
amean <- function(r,K){
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
  } else if (r < (1/(K-1)) & r != 0 & r != -1 & r != -Inf){
    a <-n(r,K)
  }
  return(a)
}


#par(mfrow=c(1,3))

#R<-seq(-10,10,by=0.1)
#K=5000000# for r from -1.5 to 1
#K=900000 # other wise
#A=c()
#B=c()
#for (i in 1:length(R)){
#  A[i]=amean(R[i],K)
#  print(i)
#}
#R1<-seq(-10,-1.1,by=0.1)
#for (i in 1:length(R1)){
#  B[i]<-(R1[i]/(R1[i]+1))*K^(1+(1/R1[i]))
#}
#B[length(R1)+1]=log(K)
#R2<-seq(-0.9,-0.1,by=0.1)
#for (i in 1:length(R2)){
# B[length(R1)+1+i]<-(R2[i]+1)^(1/R2[i])
#}
#B[length(R1)+1+length(R2)+1]=exp(1)

#R3<- seq(0.1,10,by=0.1)
#for (i in 1:length(R3)){
#  B[length(R1)+2+length(R2)+i]<-(R3[i]+1)^(1/R3[i])
#}

#plot(R,A,main="K=20000",ylab=expression(a[r]/asymptotic~result),xlab="r",type="l")
#lines(R,B, col="red",type="l")



##############################################################################################
##############################################################################################
########## The averaging method under arbitrary dependence and independence###################

#' @title The Generalized Mean Merging Function
#' @description A function to merge arbitrarily dependent p-values or independent p-values via the generalized mean method.
#' @param p numeric vector of p-values.
#' @param r the exponent of the generalized mean function in the set \eqn{[-10,\infty]\cup\{-\infty\}}; by default, r = 0.
#' @param dependence the dependence assumption of p-values, "A" or "I"; by default, dependence = "A".
#' @param censor logical; if TRUE, p-values are left-censored at 0 and right-censored at 1; by default, censor = TRUE.
#' @details  The generalized mean merging method combines p-values through the generalized mean function
#' \deqn{M_{r,K}(p_{1},\dots,p_{K})=\left(\frac{1}{K}\sum_{k=1}^{K}p_{k}^{r}\right)^{\frac{1}{r}},}
#' including the limiting cases:
#' \deqn{M_{-\infty,K}(p_{1},\dots,p_{K})=\min\{p_{1},\dots,p_{K}\};}
#' \deqn{M_{0,K}(p_{1},\dots,p_{K})=\left(\prod_{k=1}^{K}p_{k}\right)^{\frac{1}{K}};}
#' \deqn{M_{\infty,K}(p_{1},\dots,p_{K})=\max\{p_{1},\dots,p_{K}\}.}
#' @details Some special cases of the generalized mean merging method are:
#' \enumerate{
#' \item
#' The Bonferroni method for arbitrarily dependent p-values [\eqn{r=-\infty}, method="A"];
#' \item
#' The Fisher method for independent p-values [\eqn{r=0} and method="I"];
#' \item
#' The harmonic mean method for independent p-values [\eqn{r=-1}, method="I"] \insertCite{W}{pmerge};
#' \item
#' The harmonic mean method for arbitrarily dependent p-values [\eqn{r=-1}, method="A"] \insertCite{VW}{pmerge}.
#' }
#' @details  "A" dependence [dependence="A"]: The method is valid for dependent p-values. The merged p-value is calculated by
#' \deqn{F_{r,K}(p_{1},\dots,p_{K})=b_{r,K}M_{r,K}(p_{1},\dots,p_{K}),}
#' where \eqn{M_{r,K}} is the generalized mean of the p-values with exponent \eqn{r}, and \eqn{b_{r,K}} is a constant; the constant \eqn{b_{r,K}} can be found in \insertCite{VW;textual}{pmerge} and \insertCite{VWW;textual}{pmerge}.
#' @details  "I" dependence [dependence="I"]: The method is valid for independent p-values. The merged p-value is approximated by the classic central limit theorem or the generalized central limit theorem; see \insertCite{C;textual}{pmerge} for more details.  Also see \strong{Note} for the use of this method.
#' @details  If \eqn{r\le 0}, any zero p-value will produce a zero merged p-value.
#' @details  A warning is given if any of the p-values is less than zero or greater than one.
#' @import stabledist
#' @import stats
#' @importFrom Rdpack reprompt
#' @return The function returns the merged p-value.
#' @note
#' The method for independent p-values [dependence="I"]: The (generalized) central limit theorem has generally good approximations for merging independent p-values except for \eqn{r} close to \eqn{-0.5}. Therefore, we do not recommend to use \eqn{r\in(-1,0)} for independent p-values. For other choices of \eqn{r}, the number of p-values is suggested to be greater than \eqn{100}.
#' @references
#' \insertRef{C}{pmerge}
#'
#' \insertRef{VW}{pmerge}
#'
#' \insertRef{VWW}{pmerge}
#'
#' \insertRef{W}{pmerge}
#' @export

pmean <- function( p, r = 0, dependence = "A", censor = TRUE){
  K=length(p)

  if (censor == TRUE){
    p[p<0]=0
    p[p>1]=1
  } else {
    p = p}

  if (min(p) < 0 | max(p) > 1){
    warning("At least one of the input p-values are less than zero or greater than one")
  }

  if (r < -10 & r > (- Inf)){
    stop("r should be greater than -10 or equal to -Inf")
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
      pp = pnorm((gmean(p,r)^(r)-mu)*(sqrt(K/sig)))
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

      if (r >= (-1/2)){
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

#pmean(p, r=-2, dependence = "A",censor=TRUE)
#pmean(p, r=-2, dependence = "I",censor=TRUE)
