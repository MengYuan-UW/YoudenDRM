#'YoudenDRM
#'
#' @description  The Youden index is one of popular summary statistics of the receiver operating characteristic curve.
#' It has the advantage of providing a criterion to choose the “optimal” cutoff.
#' \insertCite{yuan2021semiparametric;textual}{YoudenDRM} proposed a semiparametric method to enhance the statistcal inference on the Younde index and the optimal cutoff point by using the density ratio model (DRM).
#' The proposed method covers both cases with and without a lower limit of detection.
#' @description The YoudenDRM-package contains functions to make inference on the Younde index and the optimal cutoff point under the DRM, including point estimation and confidence interval construction.
#'
#' @details
#' Let \eqn{F_0} and \eqn{F_1} denote the cumulative distribution functions of the healthy population and the diseased population, respectively.
#' The Youden index is defined as
#' \deqn{J = \max_x{F_0(x) - F_1(x)} = F_0(c) - F_1(c),}
#' where \eqn{c} is the corresponding cutoff point.
#'
#' To incorporate the information from both samples,  \insertCite{yuan2021semiparametric}{YoudenDRM} porposed to link \eqn{F_0} and \eqn{F_1} by a DRM:
#' \deqn{dF_1(x) = \exp\{\alpha + \boldsymbol{\beta}^\top \boldsymbol{q}(x)\}dF_0(x),}
#' where \eqn{dF_i} denote the density of \eqn{F_i}; \eqn{\alpha} and \eqn{\boldsymbol{\beta}} are unknown parameters for the DRM;
#' \eqn{\boldsymbol{q}(x)} is a pre-specified, non-trivial function of dimension \eqn{d};
#'   the baseline distribution \eqn{F_0} is unspecified.
#'
#'@details
#' This package contains the following material:
#' \itemize{
#' \item The Density Ratio Model (\code{\link{DRM}}): a function to fit the DRM.
#' \item The inference on the Youden index and cutoff point (\code{\link{Youden}}): a function to estimate the Youden index and cutoff point as well as construct their confidence intervals.
#' \item The goodnees-of-fit test for the DRM (\code{\link{goodnessFit}}): a test to check the validity of the DRM with a pre-specified basis function.
#' }

#' @references
#'
#' \insertRef{yuan2021semiparametric}{YoudenDRM}
#'
#' @docType package
#' @name YoudenDRM-package
NULL
