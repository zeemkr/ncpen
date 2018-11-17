#' ncpen: A package for non-convex penalized estimation for generalized linear models
#'
#' This package fits the generalized linear models with various non-convex penalties.
#' Supported regression models are Gaussian (linear), binomial Logit (logistic), multinomial Logit,
#' Poisson and Cox proportional hazard.
#' A unified algorithm is implemented in \bold{ncpen} based on the convex concave procedure
#' or difference convex algorithm that can be applied to most of existing non-convex penalties.
#' The available penalties in the package are
#' the least absolute shrinkage and selection operator(LASSO),
#' smoothly clipped absolute deviation (SCAD),
#' minimax concave penalty (MCP),
#' truncated \eqn{\ell_1}-penalty (TLP),
#' clipped LASSO (CLASSO),
#' sparse bridge (SRIDGE),
#' modified bridge (MBRIDGE),
#' and modified log (MLOG) penalties.
#'
#' The package accepts a design matrix \eqn{X} and vector of responses \eqn{y},
#' and produces the regularization path over a grid of values for the tuning parameter \code{lambda}.
#' Also provides user-friendly processes for plotting, selecting tuning parameters using cross-validation or generalized information criterion (GIC),
#' \eqn{\ell_2}-regularization, penalty weights, standardization and intercept.
#'
#' @docType package
#' @name ncpen-package
#' @note
#' This research is funded by Julian Virtue Professorship from Center for Applied Research at Pepperdine
#' Graziadio Business School and the National Research Foundation of Korea (NRF) funded
#' by Korean government (No. 2017R1C1B2010113 and 2017R1D1A1B03031239).
#'
#' @author
#' Dongshin Kim, Sunghoon Kwon and Sangin Lee
#'
#' @references
#'
#' Kim, D., Lee, S. and Kwon, S. (2018). A unified algorithm for the non-convex penalized estimation: The \code{ncpen} package.
#' \emph{http://arxiv.org/abs/1811.05061}.
#'
#' Kwon, S., Lee, S. and Kim, Y. (2016). Moderately clipped LASSO. \emph{Computational Statistics and Data Analysis}, \bold{92C}, 53-67.
#'
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems. \emph{Computational Statistics and Data Analysis}, \bold{94}, 275-286.
#'
#' Choi, H., Kim, Y. and Kwon, S. (2013). Sparse bridge estimation with a diverging number of parameters. \emph{Statistics and Its Interface}, \bold{6}, 231-242.
#'
#'
#'
NULL
