
#' @title
#' Automatic model selection regression.
#'
#'
#' @description
#' Performs k-fold cross-validation for nonconvex penalized regression models
#' to produce a regession model.
#' This is a wapper function of \code{\link{cv.ncpen}}. For advanced options, use \code{\link{cv.ncpen}}.
#'
#' @param y.vec (numeric vector) response vector.
#' @param x.mat (numeric matrix) design matrix. Each row is an observation vector.
#' @param family (character) regression model. Supported models are \code{linear} (\code{gaussian}),
#' \code{logit} (\code{binomial}) and \code{poisson}. Default is \code{linear}.
#' @param penalty (character) penalty function. Supported penalties are
#' \code{scad} (smoothly clipped absolute deviation),
#' \code{mcp} (minimax concave penalty),
#' \code{tlp} (truncated LASSO penalty),
#' \code{lasso} (least absolute shrinkage and selection operator),
#' \code{classo} (clipped LASSO),
#' \code{sridge} (sparse ridge),
#' \code{mbridge} (modified bridge) and
#' \code{mlog} (modified log).
#' Default is \code{scad}.
#' @param n.fold (numeric) the number of folds. Default value is 10. It should be 3 or greater.
#' @param intercept (logical) whether to include an intercept in the model. Default value is \code{TRUE}.
#' @param ... other parameters are same as in \code{\link{ncpen}}.
#'
#'
#' @details
#' The function runs the \code{ncpen} function for \code{n.fold+1} times.
#' The first run is to get the sequence of \code{lambda} and then the rest runs are to compute the fit with each of the folds omitted. It provides the cross validated-error based on the squared-error loss and the deviance loss.
#'
#'
#' @return An object with S3 class \code{cv.ncpen}.
#'   \item{ncpen.fit}{the fitted \code{ncpen} object.}
#'   \item{opt.ebeta}{the optimal coefficients vector selected by using the squared-error loss in the cross-validation.}
#'   \item{opt.dbeta}{the optimal coefficients vector selected by using the deviance loss in the cross-validation.}
#'   \item{cv.error}{the averaged cross-validated error for each value of \code{lambda}s.}
#'   \item{cv.deviance}{the averaged cross-validated deviance for each value of \code{lambda}s.}
#'   \item{elambda}{the \code{lambda} sequence used for computing cv error.}
#'   \item{dlambda}{the \code{lambda} sequence used for computing cv deviance.}
#'   \item{opt.elambda}{the optimal value of \code{lambda} based on cv error.}
#'   \item{opt.dlambda}{the optimal value of \code{lambda} based on cv deviance.}
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#'
#'
#' @references
#' Kwon, S., Lee, S. and Kim, Y. (2016). Moderately clipped LASSO. \emph{Computational Statistics and Data Analysis}, \bold{92C}, 53-67.
#'
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems. \emph{Computational Statistics and Data Analysis}, \bold{94}, 275-286.
#'
#' Choi, H., Kim, Y. and Kwon, S. (2013). Sparse bridge estimation with a diverging number of parameters. \emph{Statistics and Its Interface}, \bold{6}, 231-242.
#'
#'
#' @seealso
#' \code{\link{cv.ncpen}}, \code{\link{plot.cv.ncpen}}, \code{\link{coef.cv.ncpen}}, \code{\link{ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' cvfit = auto.reg(y.vec=y.vec,x.mat=x.mat,family="linear",n.fold=10)
#' coef(cvfit)
#' plot(cvfit)
#'
#' fit = cvfit$ncpen.fit
#' opt = which(cvfit$opt.elambda==fit$lambda)
#' coef(fit)[,opt] # same as coef(cvfit)
#'
#' @export
auto.reg = function(y.vec, x.mat,
                    family=c("linear","logit","poisson"),
                    penalty=c("scad","mcp","tlp","lasso","classo","sridge","mbridge","mlog"),
                    n.fold=10, df.max = 50, intercept=TRUE) {

     fam = tolower(family);
     pen = tolower(penalty);

     if(fam == "linear") {
          fam = "gaussian";
     } else if(fam == "logit") {
          fam = "binomial";
     }

     if(!is.matrix(x.mat)) {
          x.mat = as.matrix(x.mat);
     }

     cvncp = cv.ncpen(y.vec = y.vec, x.mat = x.mat, family = fam, penalty = pen, n.fold = n.fold, df.max = df.max, intercept = intercept);

     return (cvncp);
}

