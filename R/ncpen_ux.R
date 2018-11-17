#' @title
#' ncpen.reg: nonconvex penalized estimation
#' @description
#' Fits generalized linear models by penalized maximum likelihood estimation.
#' The coefficients path is computed for the regression model over a grid of the regularization parameter \code{lambda}.
#' Fits Gaussian (linear), binomial Logit (logistic), Poisson, multinomial Logit regression models, and
#' Cox proportional hazard model with various non-convex penalties.
#' @param formula (formula) regression formula. To include/exclude intercept, use \code{intercept} option
#' instead of using the "0 +" option in the formula.
#' The y value must be 0,1 for \code{binomial} and 1,2,..., for \code{multinomial}.
#' @param data (numeric matrix or data.frame) contains both y and X. Each row is an observation vector.
#' The censoring indicator must be included at the last column of the data for \code{cox}.
#' @param family (character) regression model. Supported models are
#' \code{gaussian} (or \code{linear}),
#' \code{binomial} (or \code{logit}),
#' \code{poisson},
#' \code{multinomial},
#' and \code{cox}.
#' Default is \code{gaussian}.
#' @param penalty (character) penalty function.
#' Supported penalties are
#' \code{scad} (smoothly clipped absolute deviation),
#' \code{mcp} (minimax concave penalty),
#' \code{tlp} (truncated LASSO penalty),
#' \code{lasso} (least absolute shrinkage and selection operator),
#' \code{classo} (clipped lasso = mcp + lasso),
#' \code{ridge} (ridge),
#' \code{sridge} (sparse ridge = mcp + ridge),
#' \code{mbridge} (modified bridge) and
#' \code{mlog} (modified log).
#' Default is \code{scad}.
#' @param x.standardize (logical) whether to standardize \code{x.mat} prior to fitting the model (see details).
#' The estimated coefficients are always restored to the original scale.
#' @param intercept (logical) whether to include an intercept in the model.
#' @param lambda (numeric vector) user-specified sequence of \code{lambda} values.
#' Default is supplied automatically from samples.
#' @param n.lambda (numeric) the number of \code{lambda} values.
#' Default is 100.
#' @param r.lambda (numeric) ratio of the smallest \code{lambda} value to largest.
#' Default is 0.001 when n>p, and 0.01 for other cases.
#' @param w.lambda (numeric vector) penalty weights for each coefficient (see references).
#' If a penalty weight is set to 0, the corresponding coefficient is always nonzero.
#' @param gamma (numeric) additional tuning parameter for controlling shrinkage effect of \code{classo} and \code{sridge} (see references).
#' Default is half of the smallest \code{lambda}.
#' @param tau (numeric) concavity parameter of the penalties (see reference).
#' Default is 3.7 for \code{scad}, 2.1 for \code{mcp}, \code{classo} and \code{sridge}, 0.001 for \code{tlp}, \code{mbridge} and \code{mlog}.
#' @param alpha (numeric) ridge effect (weight between the penalty and ridge penalty) (see details).
#' Default value is 1. If penalty is \code{ridge} and \code{sridge} then \code{alpha} is set to 0.
#' @param df.max (numeric) the maximum number of nonzero coefficients.
#' @param cf.max (numeric) the maximum of absolute value of nonzero coefficients.
#' @param proj.min (numeric) the projection cycle inside CD algorithm (largely internal use. See details).
#' @param add.max (numeric) the maximum number of variables added in CCCP iterations (largely internal use. See references).
#' @param niter.max (numeric) maximum number of iterations in CCCP.
#' @param qiter.max (numeric) maximum number of quadratic approximations in each CCCP iteration.
#' @param aiter.max (numeric) maximum number of iterations in CD algorithm.
#' @param b.eps (numeric) convergence threshold for coefficients vector.
#' @param k.eps (numeric) convergence threshold for KKT conditions.
#' @param c.eps (numeric) convergence threshold for KKT conditions (largely internal use).
#' @param cut (logical) convergence threshold for KKT conditions  (largely internal use).
#' @param local (logical) whether to use local initial estimator for path construction. It may take a long time.
#' @param local.initial (numeric vector) initial estimator for \code{local=TRUE}.
#' @details
#' The sequence of models indexed by \code{lambda} is fit
#' by using concave convex procedure (CCCP) and coordinate descent (CD) algorithm (see references).
#' The objective function is \deqn{ (sum of squared residuals)/2n + [alpha*penalty + (1-alpha)*ridge] } for \code{gaussian}
#' and \deqn{ (log-likelihood)/n - [alpha*penalty + (1-alpha)*ridge] } for the others,
#' assuming the canonical link.
#' The algorithm applies the warm start strategy (see references) and tries projections
#' after \code{proj.min} iterations in CD algorithm, which makes the algorithm fast and stable.
#' \code{x.standardize} makes each column of \code{x.mat} to have the same Euclidean length
#' but the coefficients will be re-scaled into the original.
#' In \code{multinomial} case, the coefficients are expressed in vector form. Use \code{\link{coef.ncpen}}.
#' @return An object with S3 class \code{ncpen}.
#'   \item{y.vec}{response vector.}
#'   \item{x.mat}{design matrix.}
#'   \item{family}{regression model.}
#'   \item{penalty}{penalty.}
#'   \item{x.standardize}{whether to standardize \code{x.mat=TRUE}.}
#'   \item{intercept}{whether to include the intercept.}
#'   \item{std}{scale factor for \code{x.standardize}.}
#'   \item{lambda}{sequence of \code{lambda} values.}
#'   \item{w.lambda}{penalty weights.}
#'   \item{gamma}{extra shrinkage parameter for \code{classo} and {sridge} only.}
#'   \item{alpha}{ridge effect.}
#'   \item{local}{whether to use local initial estimator.}
#'   \item{local.initial}{local initial estimator for \code{local=TRUE}.}
#'   \item{beta}{fitted coefficients. Use \code{coef.ncpen} for \code{multinomial} since the coefficients are represented as vectors.}
#'   \item{df}{the number of non-zero coefficients.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties.
#' \emph{Journal of the American statistical Association}, 96, 1348-60.
#' Zhang, C.H. (2010). Nearly unbiased variable selection under minimax concave penalty.
#' \emph{The Annals of statistics}, 38(2), 894-942.
#' Shen, X., Pan, W., Zhu, Y. and Zhou, H. (2013). On constrained and regularized high-dimensional regression.
#' \emph{Annals of the Institute of Statistical Mathematics}, 65(5), 807-832.
#' Kwon, S., Lee, S. and Kim, Y. (2016). Moderately clipped LASSO.
#' \emph{Computational Statistics and Data Analysis}, 92C, 53-67.
#' Kwon, S. Kim, Y. and Choi, H.(2013). Sparse bridge estimation with a diverging number of parameters.
#' \emph{Statistics and Its Interface}, 6, 231-242.
#' Huang, J., Horowitz, J.L. and Ma, S. (2008). Asymptotic properties of bridge estimators in sparse high-dimensional regression models.
#' \emph{The Annals of Statistics}, 36(2), 587-613.
#' Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models.
#' \emph{Annals of statistics}, 36(4), 1509.
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{coef.ncpen}}, \code{\link{plot.ncpen}}, \code{\link{gic.ncpen}}, \code{\link{predict.ncpen}}, \code{\link{cv.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=5,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="gaussian")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' data = cbind(y.vec, x.mat)
#' colnames(data) = c("y", paste("xv", 1:ncol(x.mat), sep = ""))
#' fit1 = ncpen.reg(formula = y ~ xv1 + xv2 + xv3 + xv4 + xv5, data = data,
#'                  family="gaussian", penalty="scad")
#' fit2 = ncpen(y.vec=y.vec,x.mat=x.mat);
#'
#' @export
ncpen.reg = function(formula,data,
                 family=c("gaussian","linear","binomial","logit","multinomial","cox","poisson"),
                 penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                 x.standardize=TRUE,intercept=TRUE,
                 lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                 df.max=50,cf.max=100,proj.min=10,add.max=10,niter.max=30,qiter.max=10,aiter.max=100,
                 b.eps=1e-7,k.eps=1e-4,c.eps=1e-6,cut=TRUE,local=FALSE,local.initial=NULL){

     family = match.arg(family);
     if(!is.data.frame(data)) {
          data = data.frame(data);
     }
     # compse y.vec and x.mat ---------------------------------
     cox.censor = NULL;
     if(family == "cox") {
          cox.censor = data[, ncol(data)];
          data = data[, -ncol(data)];
     }

     ncp.data = make.ncpen.data(formula, data);
     y.vec = ncp.data$y.vec;
     x.mat = ncp.data$x.mat;

     if(family == "cox") {
          x.mat = cbind(x.mat, cox.censor);
     }
     #---------------------------------------------------------

     ncp = ncpen(y.vec, x.mat, family, penalty,
                   x.standardize,intercept,
                   lambda,n.lambda,r.lambda,w.lambda,gamma,tau,alpha,
                   df.max,cf.max,proj.min,add.max,niter.max,qiter.max,aiter.max,
                   b.eps,k.eps,c.eps,cut,local,local.initial);

     ncp$formula = formula;
     return (ncp);
}

######################################################### cv.ncpen
#' @title
#' cv.ncpen: cross validation for \code{ncpen}
#' @description
#' performs k-fold cross-validation (CV) for nonconvex penalized regression models
#' over a sequence of the regularization parameter \code{lambda}.
#' @param formula (formula) regression formula. To include/exclude intercept, use \code{intercept} option
#' instead of using the "0 +" option in the formula.
#' The y value must be 0,1 for \code{binomial} and 1,2,..., for \code{multinomial}.
#' @param data (numeric matrix or data.frame) contains both y and X. Each row is an observation vector.
#' The censoring indicator must be included at the last column of the data for \code{cox}.
#' @param family (character) regression model. Supported models are
#' \code{gaussian} (or \code{linear}),
#' \code{binomial} (or \code{logit}),
#' \code{poisson},
#' \code{multinomial},
#' and \code{cox}.
#' Default is \code{gaussian}.
#' @param penalty (character) penalty function.
#' Supported penalties are
#' \code{scad} (smoothly clipped absolute deviation),
#' \code{mcp} (minimax concave penalty),
#' \code{tlp} (truncated LASSO penalty),
#' \code{lasso} (least absolute shrinkage and selection operator),
#' \code{classo} (clipped lasso = mcp + lasso),
#' \code{ridge} (ridge),
#' \code{sridge} (sparse ridge = mcp + ridge),
#' \code{mbridge} (modified bridge) and
#' \code{mlog} (modified log).
#' Default is \code{scad}.
#' @param x.standardize (logical) whether to standardize \code{x.mat} prior to fitting the model (see details).
#' The estimated coefficients are always restored to the original scale.
#' @param intercept (logical) whether to include an intercept in the model.
#' @param lambda (numeric vector) user-specified sequence of \code{lambda} values.
#' Default is supplied automatically from samples.
#' @param n.lambda (numeric) the number of \code{lambda} values.
#' Default is 100.
#' @param r.lambda (numeric) ratio of the smallest \code{lambda} value to largest.
#' Default is 0.001 when n>p, and 0.01 for other cases.
#' @param w.lambda (numeric vector) penalty weights for each coefficient (see references).
#' If a penalty weight is set to 0, the corresponding coefficient is always nonzero.
#' @param gamma (numeric) additional tuning parameter for controlling shrinkage effect of \code{classo} and \code{sridge} (see references).
#' Default is half of the smallest \code{lambda}.
#' @param tau (numeric) concavity parameter of the penalties (see reference).
#' Default is 3.7 for \code{scad}, 2.1 for \code{mcp}, \code{classo} and \code{sridge}, 0.001 for \code{tlp}, \code{mbridge} and \code{mlog}.
#' @param alpha (numeric) ridge effect (weight between the penalty and ridge penalty) (see details).
#' Default value is 1. If penalty is \code{ridge} and \code{sridge} then \code{alpha} is set to 0.
#' @param df.max (numeric) the maximum number of nonzero coefficients.
#' @param cf.max (numeric) the maximum of absolute value of nonzero coefficients.
#' @param proj.min (numeric) the projection cycle inside CD algorithm (largely internal use. See details).
#' @param add.max (numeric) the maximum number of variables added in CCCP iterations (largely internal use. See references).
#' @param niter.max (numeric) maximum number of iterations in CCCP.
#' @param qiter.max (numeric) maximum number of quadratic approximations in each CCCP iteration.
#' @param aiter.max (numeric) maximum number of iterations in CD algorithm.
#' @param b.eps (numeric) convergence threshold for coefficients vector.
#' @param k.eps (numeric) convergence threshold for KKT conditions.
#' @param c.eps (numeric) convergence threshold for KKT conditions (largely internal use).
#' @param cut (logical) convergence threshold for KKT conditions  (largely internal use).
#' @param local (logical) whether to use local initial estimator for path construction. It may take a long time.
#' @param local.initial (numeric vector) initial estimator for \code{local=TRUE}.
#' @param n.fold (numeric) number of folds for CV.
#' @param fold.id (numeric vector) fold ids from 1 to k that indicate fold configuration.
#' @details
#' Two kinds of CV errors are returned: root mean squared error and negative log likelihood.
#' The results depends on the random partition made internally.
#' To choose an optimal coefficients form the cv results, use \code{\link{coef.cv.ncpen}}.
#' \code{ncpen} does not search values of \code{gamma}, \code{tau} and \code{alpha}.
#' @return An object with S3 class \code{cv.ncpen}.
#'   \item{ncpen.fit}{ncpen object fitted from the whole samples.}
#'   \item{fold.index}{fold ids of the samples.}
#'   \item{rmse}{rood mean squared errors from CV.}
#'   \item{like}{negative log-likelihoods from CV.}
#'   \item{lambda}{sequence of \code{lambda} used for CV.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties.
#' \emph{Journal of the American statistical Association}, 96, 1348-60.
#' Zhang, C.H. (2010). Nearly unbiased variable selection under minimax concave penalty.
#' \emph{The Annals of statistics}, 38(2), 894-942.
#' Shen, X., Pan, W., Zhu, Y. and Zhou, H. (2013). On constrained and regularized high-dimensional regression.
#' \emph{Annals of the Institute of Statistical Mathematics}, 65(5), 807-832.
#' Kwon, S., Lee, S. and Kim, Y. (2016). Moderately clipped LASSO.
#' \emph{Computational Statistics and Data Analysis}, 92C, 53-67.
#' Kwon, S. Kim, Y. and Choi, H.(2013). Sparse bridge estimation with a diverging number of parameters.
#' \emph{Statistics and Its Interface}, 6, 231-242.
#' Huang, J., Horowitz, J.L. and Ma, S. (2008). Asymptotic properties of bridge estimators in sparse high-dimensional regression models.
#' \emph{The Annals of Statistics}, 36(2), 587-613.
#' Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models.
#' \emph{Annals of statistics}, 36(4), 1509.
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{plot.cv.ncpen}}, \code{\link{coef.cv.ncpen}}, \code{\link{ncpen}}, \code{\link{predict.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=5,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="gaussian")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' data = cbind(y.vec, x.mat)
#' colnames(data) = c("y", paste("xv", 1:ncol(x.mat), sep = ""))
#' fit1 = cv.ncpen.reg(formula = y ~ xv1 + xv2 + xv3 + xv4 + xv5, data = data, n.lambda=10,
#'                     family="gaussian", penalty="scad")
#' fit2 = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="gaussian", penalty="scad")
#' coef(fit1)
#'
#' @export
cv.ncpen.reg = function(formula,data,
                    family=c("gaussian","linear","binomial","logit","multinomial","cox","poisson"),
                    penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                    x.standardize=TRUE,intercept=TRUE,
                    lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                    df.max=50,cf.max=100,proj.min=10,add.max=10,niter.max=30,qiter.max=10,aiter.max=100,
                    b.eps=1e-6,k.eps=1e-4,c.eps=1e-6,cut=TRUE,local=FALSE,local.initial=NULL,
                    n.fold=10,fold.id=NULL){

     family = match.arg(family);
     if(!is.data.frame(data)) {
          data = data.frame(data);
     }
     # compse y.vec and x.mat ---------------------------------
     cox.censor = NULL;
     if(family == "cox") {
          cox.censor = data[, ncol(data)];
          data = data[, -ncol(data)];
     }

     ncp.data = make.ncpen.data(formula, data);
     y.vec = ncp.data$y.vec;
     x.mat = ncp.data$x.mat;

     if(family == "cox") {
          x.mat = cbind(x.mat, cox.censor);
     }

     #---------------------------------------------------------
     cvncp = cv.ncpen(y.vec, x.mat, family, penalty,
                   x.standardize,intercept,
                   lambda,n.lambda,r.lambda,w.lambda,gamma,tau,alpha,
                   df.max,cf.max,proj.min,add.max,niter.max,qiter.max,aiter.max,
                   b.eps,k.eps,c.eps,cut,local,local.initial,
                   n.fold,fold.id);

     cvncp$formula = formula;
     return (cvncp);
}

