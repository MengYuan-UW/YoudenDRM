DRMest = function(x,y,qt,r = NULL,totalSize, method = "optimal"){
  model = DRM(x,y,qt,r,totalSize,method)
  coef = model$coef
  t = model$sample_obs
  zeta = model$prob_obs
  P0 = model$P0
  P1 = model$P1
  #ecdf
  F0_hat = function(c){1 - zeta[1] + sum(P0*ifelse(t <= c,1,0))}
  F1_hat = function(c){1 - zeta[2] + sum(P1*ifelse(t <= c,1,0))}

  if(is.null(r)){
    value  = t
    }else{
      value = c(r,t)
    }
  F0 = sapply(sort(value),F0_hat)
  F1 = sapply(sort(value),F1_hat)
  #reorder = order(t)
  #F00 = (1- zeta[1]) + cumsum(P0[reorder])
  #F1 = (1- zeta[2])+ cumsum(P1[reorder])
  dif = F0 - F1
  #calulate J and c
  f = function(t){eval(parse(text = paste(coef,c(1,qt),sep = "*",collapse = "+")))}
  points = tryCatch(uniroot.all(f,lower = min(value),upper = max(value)),
                    error = function(e) NULL)
  if(is.null(points)){
    J = max(dif)
    cut = sort(value)[dif == J]
    c = mean(cut)
    warning("The cutoff point is out of range")
  }else{
    index = findInterval(points, sort(value))
    dif_J= dif[index]
    J = max(dif_J)
    cut = points[index == index[dif_J == J]]
    c = mean(cut)
  }
  #preperation for variance estimate
  qt_data = lapply(1:length(qt),function(i){eval(parse(text = qt[i]))})
  qt_matrix = cbind(1,matrix(unlist(qt_data),nrow = length(t),byrow = F))

  if(is.null(r)){
    n = length(y)+length(x)
    rho = length(y)/length(x)
  }else{
    if(missing(totalSize)){
      stop("Please provide the total sample sizes")
    }else{
      n=sum(totalSize)
      rho =totalSize[2]/totalSize[2]
    }
  }

  delta = 1/(1+rho*exp(qt_matrix%*%coef))

  qt1_matrix = qt_matrix[,-1]
  qt2_list = lapply(1:length(qt),function(i){
    phase = paste(qt[i],qt,sep = "*",collapse = ",")
    eval(parse(text = paste("cbind(",phase,")")))
  })
  qt2_matrix = do.call("cbind",qt2_list)

  qc = gsub("t","c",qt)
  expres = parse(text = paste(coef[-1],qc,sep = "*",collapse = "+"))
  dqc = eval(D(expres,"c"))
  qc = c(1,unlist(lapply(1:length(qc),
                         function(i){eval(parse(text = qc[i]))})))
  #### for sigma_c
  A0 = sum(delta*P1)
  A1 = t(qt1_matrix)%*%(delta*P1)
  A2 = matrix(t(qt2_matrix)%*%(delta*P1),nrow = length(A1),byrow = T)
  A = rbind(c(A0,A1),cbind(A1,A2))
  S = rho/(1+rho)*A
  V = S- rho*c(A0,A1)%*%t(c(A0,t(A1)))
  S0 = solve(S)%*%V%*%solve(S)
  Sigma2_c = t(qc)%*%S0%*%qc/(dqc)^2
  se_c = sqrt(Sigma2_c/n)

  ######## for J ###########
  F0_c = 1 - zeta[1] + sum(P0*ifelse(t <= c,1,0))
  F1_c = 1 - zeta[2] + sum(P1*ifelse(t <= c,1,0))
  A0_c = sum(delta*P1*ifelse(t <= c,1,0))
  A1_c = t(qt1_matrix)%*%(delta*P1*ifelse(t <= c,1,0))
  Sigma2_J = (rho+1)*F0_c*(1-F0_c)+(rho+1)/rho*F1_c*(1-F1_c) -
    (rho+1)^3/rho*(A0_c - t(c(A0_c,A1_c))%*%solve(A)%*%c(A0_c,A1_c))
  se_J = sqrt(Sigma2_J/n)
  list(J = J, sdJ = as.numeric(se_J), c = c, sdc = as.numeric(se_c))
}

##################################################

#' @title Inference on the Youden index and optimal cutoff point
#' @description  Estimate the Youden index and its corrreponding cutoff point under the density ratio model (DRM) as well as construct confidence intervals based on the provided data.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param qt pre-sepecified basis functions of t in the exponential term. For example, \code{qt = c("t","log(t)")}.
#' @param r  the value of the lower limit of detetion (LLOD). The default case is where no LLOD exists.
#' @param totalSize the total sample sizes \eqn{(n_0,n_1)}. When there is no LLOD (\code{r = NULL}), this argument is optional. See \code{help(DRM)} `Details'.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for details.
#' @param CItype the method to be used for confidence interval construction. See `Details'.
#' @param level confidence level of the confidence interval. The default value is 0.95.
#'
#'@details Let \eqn{F_0} and \eqn{F_1} denote the cumulative distribution functions of the healthy population and the diseased population, respectively.
#' The Youden index is defined as
#' \deqn{J = \max_x{F_0(x) - F_1(x)} = F_0(c) - F_1(c),}
#' where \eqn{c} is the corresponding cutoff point.
#'
#' Under the DRM, the estimator \eqn{\hat c} of cutoff point is the solution to the equation
#' \deqn{\hat\alpha + \hat{\boldsymbol{\beta}}^\top \boldsymbol{q}(x) = 0,}
#' where \eqn{\hat\alpha} and \eqn{\hat{\boldsymbol{\beta}}} are the estimators of the paramters of DRM.
#' If multiple solution exist in the range of the observed samples, we choose the one that attains the maximum of \eqn{\hat F_0(\hat c) - \hat F_1(\hat c)} with
#' \eqn{\hat F_0} and \eqn{\hat F_1} being  the estimator of \eqn{F_0} under the DRM.
#' The explicit form of \eqn{\hat F_0} can be found in \insertCite{yuan2021semiparametric;textual}{YoudenDRM}.
#' If no solution exists  in the range of the observed samples, we find the \eqn{\hat c} by the definition.
#' With \eqn{\hat c}, the estimator of Youden index is \eqn{\hat J = \hat F_0(\hat c) - \hat F_1(\hat c)}.
#'
#'@details The argument \code{CItype} refers to difference confidence intervals of the Younde index and cutoff point
#'under the DRMs in \insertCite{yuan2021semiparametric;textual}{YoudenDRM}.
#'\itemize{
#'\item  \code{"None"}: no confidence intervals for the Younde index and cutoff point is constructed;
#'\item  \code{"NA-DRM"}: the Wald-type confidence intervals for the Younde index and cutoff point based on the normal approximation;
#'\item  \code{"logit-DRM"}: the Wald-type confidence interval for the Youden index using logit transformation.
#'The Wald-type confidence interval for the cutoff point is contructed based on the normal approximation.
#'}
#'
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item  Youden: the estimate, asymptotic standard deviation (ASD), and the confidence intervals (lower bound and uppper bound) of the Youden index.
#' When \code{CItype = "None"}, it only returns the estimate and ASD.
#' \item  cutoff: the estimate, ASD, and the confidence intervals (lower bound and uppper bound) of the optimal cutoff point.
#' When \code{CItype = "None"}, it only returns the estimate and ASD.
#' }
#'
#' @import stats
#' @import nnet
#' @import rootSolve
#'
#'
#' @examples
#' #Example 1 (without LLOD)
#' #DMD is the dataset of data application in the paper.
#' x = DMD$CK[DMD$Status == 0]
#' y = DMD$CK[DMD$Status == 1]
#' Youden(x,y,qt = "t",CItype ="logit-DRM")
#' @examples
#' #Example 2 (with a LLOD)
#' #Date generation
#' set.seed(123456)
#' x = rlnorm(50, meanlog = 2.5, sdlog = sqrt(0.09))
#' y = rlnorm(50, meanlog = 2.87,sdlog = sqrt(0.25))
#' #Create the LLOD
#' r = qlnorm(0.15,meanlog = 2.5, sdlog = sqrt(0.09))
#' #Observed samples
#' x = x[x>=r]
#' y = y[y>=r]
#' # The basis function
#' qt = c("log(t)","log(t)^2")
#' Youden(x,y,qt,r,totalSize = c(50,50),CItype ="logit-DRM")
#' @references
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#'
#' @export
Youden = function(x,y,qt,r = NULL, totalSize,method = "optimal", CItype, level = 0.95){
  model = DRMest(x,y,qt,r,totalSize,method)
  quant = qnorm((1+level)/2)
  if(missing(CItype)){CItype == "None"}
  if(CItype == "None"){
    J = c(model$J,model$sdJ)
    c = c(model$c,model$sdc)
    names(J) = c("estimate","ASD")
    names(c) = c("estimate","ASD")
  }
  if(CItype == "NA-DRM"){
    CI_c = model$c +c(-1,1)*quant*model$sdc
    CI_J= model$J +c(-1,1)*quant*model$sdJ
    J = c(model$J,model$sdJ,CI_J)
    c = c(model$c,model$sdc,CI_c)
    names(J) = c("estimate","ASD","lower bound","upper bound")
    names(c) = c("estimate","ASD","lower bound","upper bound")
  }
  if(CItype == "logit-DRM"){
    CI_c = model$c +c(-1,1)*quant*model$sdc
    CI =log(model$J) - log(1 - model$J) + c(-1,1)*quant*model$sdJ/(model$J*(1-model$J))
    CI_J =exp(CI)/(1 + exp(CI))
    J = c(model$J,model$sdJ,CI_J)
    c = c(model$c,model$sdc,CI_c)
    names(J) = c("estimate","ASD","lower bound","upper bound")
    names(c) = c("estimate","ASD","lower bound","upper bound")
  }
  return(list(Youden = J, cutoff = c))
}
