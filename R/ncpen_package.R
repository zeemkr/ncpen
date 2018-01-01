#' ncpen: A package for non-convex penalized estimations.
#'
#' \bold{ncpen} estimaties generalized linear models with various non-convex penalties.
#'
#' \bold{ncpen} estimaties generalized linear models including gaussian (linear), poisson and binomial.
#' The implemented penalties are
#' smoothly clipped absolute deviation (SCAD),
#' minimax concave penalty (MCP),
#' truncated LASSO penalty (TLP),
#' ???s clipped LASSO (CLASSO), ???e
#' ???s bridge penalty (SRIDGE), ???e
#' modified bridge penalty (MBRIDGE),
#' ???s (MLOG), ???e
#' and least absolute shrinkage and selection operator (LASSO).
#'
#' The algorithm implements recent technical advances for the optimization of penalized estimations:
#' convex-concave procedure for nonconvex optimization,
#' local quadratic approximation (LQA),
#' a descent directional modification of LQA for GLM,
#' coordinate-wise descent algorithm for quadratic \eqn{l_1}-penalized problems,
#' active set optimization for speed up of the algorithm,
#' warm-start strategy for pathwise optimization and more.
#'
#' The package also provides various functions and options for user-specific choices of an initial values,
#' penalty terms, \eqn{l_2}-regularization, observation and penalty weights, standardization and intercept.
#'
#' @docType package
#' @name ncpen-package
#' @note
#' This project is funded by Julian Virtue Professorship from
#' Center for Applied Research at Graziadio School of Business and Management at
#' Pepperdine University.
#'
#' @author
#' Dongshin Kim, Sunghoon Kwon and Sangin Lee
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#'
NULL
