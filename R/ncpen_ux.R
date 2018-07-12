
#' @title
#' Penalized Estimation Regression
#'
#'
#' @description
#' Selects a regression model by the penalized estimation.
#' This fuction runs k-fold cross-validation for nonconvex penalized regression models
#' to produce a regession model.
#' This is a wapper function of \code{\link{cv.ncpen}}. For advanced options, use \code{\link{cv.ncpen}} instead.
#'
#' @param y.vec (numeric vector) response vector.
#' @param x.mat (numeric matrix) design matrix. Each row is an observation vector.
#' @param family (character) regression model. Supported models are \code{linear} (\code{gaussian}),
#' \code{logit} (\code{binomial}), \code{multinomial}, \code{cox} and \code{poisson}. Default is \code{linear}.
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
#' @param df.max (numeric) the expected maximum number of non-zero coefficients. Default value is 50.
#' @param intercept (logical) whether to include an intercept in the model. Default value is \code{TRUE}.
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
#' #cvfit = pen.reg(y.vec=y.vec,x.mat=x.mat,family="linear",n.fold=10)
#' #coef(cvfit)
#' #plot(cvfit)
#'
#' #fit = cvfit$ncpen.fit
#' #opt = which(cvfit$opt.elambda==fit$lambda)
#' #coef(fit)[,opt] # same as coef(cvfit)
#'
#' @export
pen.reg = function(y.vec, x.mat,
                   family=c("linear","logit", "multinomial", "cox", "poisson"),
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

#' @title
#' Penalized Estimation Regression
#'
#'
#' @description
#' Selects a regression model by the penalized estimation.
#' This function performs k-fold cross-validation for nonconvex penalized regression models
#' to produce a regession model.
#' This is a wapper function of \code{\link{cv.ncpen}}. For advanced options, use \code{\link{cv.ncpen}} instead.
#'
#' @param formula (formula) regression formula. To include/exclude intercept, use \code{intercept} option
#' instead of using the "0 +" option in the formula.
#' @param data (numeric matrix or data.frame) contains both y and X. Each row is an observation vector.
#' @param family (character) regression model. Supported models are \code{linear} (\code{gaussian}),
#' \code{logit} (\code{binomial}), \code{multinomial}, \code{cox} and \code{poisson}. Default is \code{linear}.
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
#' @param df.max (numeric) the expected maximum number of non-zero coefficients. Default value is 50.
#' @param intercept (logical) whether to include an intercept in the model. Default value is \code{TRUE}.
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
#' data = data.frame(y = y.vec, x.mat)
#'
#' #cvfit = pen.reg2(formula = y ~ ., data = data, family="linear",n.fold=10, intercept = TRUE)
#' #coef(cvfit)
#' #plot(cvfit)
#'
#' #fit = cvfit$ncpen.fit
#' #opt = which(cvfit$opt.elambda==fit$lambda)
#' #coef(fit)[,opt] # same as coef(cvfit)
#'
#' @export
pen.reg2 = function(formula, data,
                    family=c("linear","logit","multinomial","cox","poisson"),
                    penalty=c("scad","mcp","tlp","lasso","classo","sridge","mbridge","mlog"),
                    n.fold=10, df.max = 50, intercept=TRUE) {

     # compse y.vec and x.mat ---------------------------------
     ncp.data = make.ncpen.data(formula, data);
     y.vec = ncp.data$y.vec;
     x.mat = ncp.data$x.mat;
     #---------------------------------------------------------

     cvncp = pen.reg(y.vec, x.mat, family, penalty, n.fold, df.max, intercept);
     cvncp$formula = formula;
     return (cvncp);
}

