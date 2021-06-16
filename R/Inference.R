#' @title Point estimation of the Youden index and optimal cutoff point
#' @description Estimate the Youden index and its corrreponding cutoff point and standard deviation of estimators under DRM based on the provided data.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param qt pre-sepecified basis functions of t in the exponential term. For example, \code{qt = c("t","log(t)")}.
#' @param r  the value of the LLOD. The default case is where no LLOD exists.
#' @param totalSize the total sample sizes \eqn{(n_0,n_1)}. When there is no LLOD (\code{r = NULL}), this argument is optional. See \code{help(DRM)} `Details'.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for details.
#'
#' @import stats
#' @import nnet
#' @import rootSolve
#'
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item J:  estimate of the Youden index,
#' \item sdJ: standard devatation of the estimator of the Youden index,
#' \item c: sdcestimate of the optimal cutoff point,
#' \item sdc: standard devatation of the estimator of the optimal cutoff point.
#'}
#' @examples
#' #DMD is the dataset of data application in \insertCite{yuan2021semiparametric;textual}{YoudenDRM}.
#' x = DMD$CK[DMD$Status == 0]
#' y = DMD$CK[DMD$Status == 1]
#' DRMest(x,y,qt = "t")
#' @examples
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
#' DRMest(x,y,qt,r,totalSize = c(50,50))
#' @references
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#' @export

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

#' @title Confidence intervals for the Youden index and optimal cutoff point
#' @description Construct CIs for the Youden index and its corrreponding cutoff point under DRMs based on the provided data.
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param qt pre-sepecified basis functions of t in the exponential term. For example, \code{qt = c("t","log(t)")}.
#' @param r  the value of the LLOD. The default case is where no LLOD exists.
#' @param totalSize the total sample sizes \eqn{(n_0,n_1)}. When there is no LLOD (\code{r = NULL}), this argument is optional. See \code{help(DRM)} `Details'.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See \code{help(DRM)} for details.
#' @param logit whether the logit transformation is applied when constructing confidence interval for the Youden index.
#' @param level confidence level of the confidence interval. The default value is 0.95.
#' @import stats
#' @import nnet
#' @import rootSolve
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item  CI_J: confidence interval of the Youden index.
#' \item  CI_c :confidence interval of the optimal cutoff point.
#' }
#' @examples
#' #DMD is the dataset of data application in \insertCite{yuan2021semiparametric;textual}{YoudenDRM}.
#' x = DMD$CK[DMD$Status == 0]
#' y = DMD$CK[DMD$Status == 1]
#' DRMci(x,y,qt = "t")
#' @examples
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
#' DRMci(x,y,qt,r,totalSize = c(50,50))
#' @references
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#' @export
DRMci = function(x,y,qt,r = NULL, totalSize,method = "optimal", logit = T, level = 0.95){
  model = DRMest(x,y,qt,r,totalSize,method)
  quant = qnorm((1+level)/2)
  CI_c = model$c +c(-1,1)*quant*model$sdc
  if(logit == T){
    CI =log(model$J) - log(1 - model$J) + c(-1,1)*quant*model$sdJ/(model$J*(1-model$J))
    CI_J =exp(CI)/(1 + exp(CI))
  }else{
    CI_J= model$J +c(0,-1,1)*quant*model$sdJ
  }
  return(list(CI_J = CI_J,CI_c = CI_c ))
}
