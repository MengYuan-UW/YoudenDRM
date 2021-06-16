
KStest = function(x,y,qt,r = NULL,totalSize,method = "optimal"){
  model = DRM(x,y,qt,r,totalSize,method)
  coef = model$coef
  t = model$sample_obs
  #reorder = order(t)
  zeta = model$prob_obs
  P0 = model$P0
  P1 = model$P1
  if(is.null(r)){
    n0=length(x)
    n1=length(y)
  }else{
    if(missing(totalSize)){
      stop("Please provide the total sample sizes: totalSize = c(n0,n1)")
    }else{
      n0=totalSize[1]
      n1=totalSize[2]
    }
  }
  rho = c(n0,n1)/(n0+n1)
  F0_hat = function(c){1 - zeta[1] + sum(P0*ifelse(t <= c,1,0))}
  F0_emp = function(c){sum(ifelse(x <= c,1,0))/length(x)}
  F1_hat = function(c){1 - zeta[2] + sum(P1*ifelse(t <= c,1,0))}
  F1_emp = function(c){sum(ifelse(y <= c,1,0))/length(y)}
  Delta0 =  function(c){sqrt(length(x)+length(y))*abs(F0_hat(c) - F0_emp(c))}
  Delta1 =  function(c){sqrt(length(x)+length(y))*abs(F1_hat(c) - F1_emp(c))}
  if(is.null(r)){
    KS_stat = c(max(sapply(t,Delta0)),max(sapply(t,Delta1)))
  }else{
    KS_stat = c(max(sapply(c(r,t),Delta0)),max(sapply(c(r,t),Delta1)))
  }
  Delta_stat = mean(rho*KS_stat)
  if(is.null(r)){
    list(stat = Delta_stat,P0 = P0, P1 = P1, t = t)
  }else{
    list(stat = Delta_stat,P0 = c(1-zeta0,P0), P1 = c(1-zeta1,P1), t = c(r,t))
  }

}

#'
#' @title Goodness-of-fit test
#' @description Goodness-of-fit test for the validity of the DRM with a pre-specified basis function \code{qt}.
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param qt pre-sepecified basis functions of t in the exponential term. For example, \code{qt = c("t","log(t)")}.
#' @param r  the value of the LLOD. The default case is where no LLOD exists.
#' @param totalSize the total sample sizes \eqn{(n_0,n_1)}. When there is no LLOD (\code{r = NULL}), this argument is optional. See \code{help(DRM)} `Details'.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for details.
#' @param B  number of boostrap reptitions.
#'
#' @details
#' \insertCite{qin1997goodness;textual}{YoudenDRM} defined a  Kolmogorov-Smirnov-type test statistic as
#' \deqn{\Delta_n = \sup_{-\infty \leq x \leq \infty}   \sqrt{n}|\hat F_0(x) - \bar F_0(x)|,}
#' where \eqn{\hat F_0} refers to the estimator of \eqn{F_0} under the DRM
#' and \eqn{\bar F_0}  represents the empirical cumulative distribution function of \eqn{F_0}.
#' The explicit form of \eqn{\hat F_0} can be found in \insertCite{yuan2021semiparametric;textual}{YoudenDRM}.
#'
#' We reject the null hypothesis that the DRM is satisfied if \eqn{\Delta_n} is greater than some critical value.
#' The limiting distribution of \eqn{\Delta_n} has a complicated form.
#' \insertCite{qin1997goodness;textual}{YoudenDRM} suggested to use a Bootstrap method to find the critical value.
#'
#' @import stats
#' @import nnet
#' @import rootSolve
#' @importFrom Rdpack reprompt
#'
#' @return The function returns a list containing:
#' \itemize{
#' \item KSstat: the test statistic,
#' \item pval: p value of the test.
#' }
#' @examples
#' #DMD is the dataset of data application in paper
#' x = DMD$CK[DMD$Status == 0]
#' y = DMD$CK[DMD$Status == 1]
#' set.seed(123456)
#' goodnessFit(x,y,qt = "t",B = 1000)
#'
#' @references
#'
#' \insertRef{qin1997goodness}{YoudenDRM}
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#' @export
#'
goodnessFit = function(x,y,qt,r = NULL,totalSize, method = "optimal",B){
  test = KStest(x,y,qt,r,totalSize,method)
  t = test$t
  P0 = test$P0
  P1 = test$P1
  n0 = length(x)
  n1 = length(y)
  KS_stat = test$stat
  samplex = matrix(sample(t, B*n0, replace=T, prob=P0),
                   nrow=length(x), ncol=B)
  sampley = matrix(sample(t, B*n1, replace=T, prob=P1),
                   nrow=length(y), ncol=B)
  boot_sample <- rbind(samplex, sampley)
  boot_test = function(sam){
    samx = sam[1:n0]
    samy = sam[-(1:n0)]
    as.numeric(KStest(samx,samy,qt,r,totalSize,method)$stat)
  }
  boot_stat <- apply(boot_sample, 2, boot_test)
  pval <- mean(na.omit(boot_stat) > KS_stat)
  list(pval = pval, KSstat = KS_stat )
}
