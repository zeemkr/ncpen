#' ncpen: A package for noncovex penalized estimations.
#'
#' The ncpen package provides ....
#' and many
#'
#' @section ncpen functions:
#' The ncpen functions ...
#'
#' @docType package
#' @name ncpen-package
#' @description
#' \code{ncpen} fits generalized linear models by using penalized maximum likelihood estimation.
#' It covers Various non-convex penalties such as SCAD and MCP for linear, logistic and poisson regression models.
#'
#' @param y.vec numeric vector; samples of dependent variable
#' @param x.mat numeric matrix; samples of independent variables
#' @param family character; model type; defalut is "gaussian"
#'
#' @return
#' This returns...... If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#'
#' @note
#' Leave some notes here.
#'
#' @author
#' Sunghoon Kwon, Dongshin Kim
#'
#' @references
#' Paper 1 by Big Name
#'
#' @seealso
#' See this also....
#'
#' @examples
#' fam = "lin"
#' pen = "scad"
#'
#' a = 3 + 4;
#' a
#'
NULL
