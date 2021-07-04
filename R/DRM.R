#' @title Density Ratio Model (DRM)
#' @description Fit a semiparametric two-sample DRM based on data with or without the lower limit of detetion (LLOD).
#'
#' @param x observed sample from \eqn{F_0}.
#' @param y observed sample from \eqn{F_1}.
#' @param qt pre-sepecified basis functions of t in the exponential term. For example, \code{qt = c("t","log(t)")}.
#' @param r  the value of the LLOD. The default case is where no LLOD exists.
#' @param totalSize the total sample sizes \eqn{(n_0,n_1)}. When there is no LLOD (\code{r = NULL}), this argument is optional. See `Details'.
#' @param method the method used to fit the DRM. The options include \code{"glm"}, \code{"multinom"} and \code{"optimal"}.
#' The default setting is \code{"optimal"}. See `Details'.
#'
#' @details Denote \eqn{\{X_{01},\ldots,X_{0n_0}\}} and \eqn{\{X_{11},\ldots,X_{1n_1}\}} as two independent random samples coming from the healthy and diseased populations, respectively.
#' Let \eqn{F_0} and \eqn{F_1} denote the cumulative distribution functions of the healthy population and the diseased population.
#' \insertCite{yuan2021semiparametric}{YoudenDRM} porposed to link \eqn{F_0} and \eqn{F_1} by a DRM:
#' \deqn{dF_1(x) = \exp\{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(x)\}dF_0(x),}
#' where \eqn{dF_i} denote the density of \eqn{F_i}; \eqn{\alpha} and \eqn{\boldsymbol{\beta}} are unknown parameters for the DRM;
#' \eqn{\boldsymbol{q}(x)} is a pre-specified, non-trivial function of dimension \eqn{d};
#'  the baseline distribution \eqn{F_0} is unspecified.
#'
#' Suppose the biomarker has a LLOD, denoted as \eqn{r}. Let \eqn{m_0} and \eqn{m_1} be the numbers of observations above the LLOD \eqn{r} in the healthy and diseased groups, respectively.
#' Without loss of generality, we assume the first \eqn{m_i} observations in sample \eqn{i}, \eqn{X_{i1},\ldots,X_{im_i}} are above the LLOD.
#' According to \insertCite{yuan2021semiparametric;textual}{YoudenDRM},
#'the estimators of \eqn{\alpha} and \eqn{\boldsymbol{\beta}} maximize the following empirical log-likelihood function
#' \deqn{\ell_{n}(\alpha,\boldsymbol{\beta} ) = \sum_{j = 1}^{m_1} \{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(X_{1j})\} - \sum_{i = 0}^1\sum_{j = 1}^{m_i}\log\left[1+\frac{n_1}{n_0} \exp\{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(X_{ij})\}\right].}
#'
#' @details Since the DRM is equivalent to the logistic regression model with some justification of \eqn{\alpha},
#' we can use function \code{glm} or \code{multinom} to solve for parameters \eqn{\alpha} and \eqn{\boldsymbol{\beta}} instead of maximizing the above likelihood function.
#' The method \code{optimal} estimates the parameters \eqn{\alpha} and \eqn{\boldsymbol{\beta}} by \code{glm} or \code{multinom}, whichever gives the larger likelihood \eqn{\ell_{n}(\alpha,\boldsymbol{\beta})}.
#'
#' @return The function returns a list containing the following components:
#' \itemize{
#' \item prob_obs: the estimates of \eqn{\zeta_0 = P(x_{01}\geq r)} and \eqn{\zeta_1 = P(x_{11}\geq r)}. If there is no LLOD, it returns (1,1).
#' \item sample_obs: the combined observed samples.
#' \item coef: the estimates of parameters \eqn{\alpha} and \eqn{\boldsymbol{\beta}}.
#' \item P0: the estimated density of \eqn{F_0} at each observed data point.
#' \item P1: the estimated density of \eqn{F_1} at each observed data point.
#' }
#' @import nnet
#' @import stats
#' @import rootSolve
#' @importFrom Rdpack reprompt
#' @references
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#'
#' @export
DRM = function(x,y,qt,r = NULL,totalSize,method = "optimal"){
  # total size of sample from f_0 and f_1
  x_obs = x
  y_obs = y

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
  # the size of sample above the r from f_0 and f_1
  m0 = length(x_obs)
  m1 = length(y_obs)

  zeta0 = m0/n0
  zeta1 = m1/n1

  group=c(rep(0,m0),rep(1,m1))
  t = c(x_obs,y_obs)
  td = sort(t)
  qt_data = lapply(1:length(qt),function(i){eval(parse(text = qt[i]))})
  qt_matrix = cbind(1,matrix(unlist(qt_data),nrow = length(t),byrow = F))
  qy = gsub("t","y",qt)
  qy_data = lapply(1:length(qy),function(i){eval(parse(text = qy[i]))})
  qy_matrix = cbind(1, matrix(unlist(qy_data),nrow = length(y),byrow = F))

  dual = function(coef){
    pi = 1/(1+n1/n0*exp(qt_matrix%*%coef))
    loglik = sum(log(pi)) + sum(qy_matrix%*%coef)
    loglik
  }

  glm_frame =data.frame(cbind(group,qt_matrix[,-1]))
  formula = as.formula(paste("group~",
                             paste(names(glm_frame)[-1],collapse = "+"),collapse = ""))
  if(method == "optimal"){
    out1 = tryCatch(multinom(formula,data = glm_frame, trace = F),
                    error = function(e) NULL)
    out2 = tryCatch(glm(formula,data = glm_frame,family = binomial(link = "logit")),
                    error = function(e) NULL)
    if (is.null(out1) & is.null(out2) ) {stop("parameters of the density ratio model can not been solved")}
    if (is.null(out1)){
      dual1 = -1e5
      outcoef1 = NULL
    }else{
      outcoef1=coef(out1)
      outcoef1[1]=outcoef1[1]-log(n1/n0)
      dual1 = dual(outcoef1)
    }
    if (is.null(out2)){
      dual2 = -1e5
      outcoef2 = NULL
    }else{
      outcoef2 = coef(out2)
      outcoef2[1]=outcoef2[1]-log(n1/n0)
      dual2 = dual(outcoef2)
    }
    # the estimated coefficients
    if (dual1 > dual2){
      outcoef = outcoef1
    }else{
      outcoef = outcoef2
    }
  }
  if(method == "glm"){
    out = tryCatch(glm(formula,data = glm_frame,family = binomial(link = "logit")),
                   error = function(e) NULL)
    if(is.null(out)){
      stop("parameters of the density ratio model can not been solved")
    }else{
      outcoef=coef(out)
      outcoef[1]=outcoef[1]-log(n1/n0)
    }
  }
  if(method == "multinom"){
    out = tryCatch(multinom(formula,data = glm_frame, trace = F),
                   error = function(e) NULL)
    if(is.null(out)){
      stop("parameters of the density ratio model can not been solved")
    }else{
      outcoef=coef(out)
      outcoef[1]=outcoef[1]-log(n1/n0)
    }
  }

  P0 = 1/(n0+n1*exp(qt_matrix%*%outcoef))
  P1 = exp(qt_matrix%*%outcoef)*P0
  list(prob_obs = c(zeta0,zeta1), sample_obs = t,
       coef = outcoef,P0 = P0,P1 = P1)
}
