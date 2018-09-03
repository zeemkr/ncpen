#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline lines plot par
#' @importFrom stats formula median model.matrix rbinom rnorm rpois runif as.formula
#' @useDynLib ncpen, .registration = TRUE
######################################################### ncpen
#' @title
#' ncpen: nonconvex penalized estimation
#' @description
#' Fits generalized linear models by penalized maximum likelihood estimation.
#' The coefficients path is computed for the regression model over a grid of the regularization parameter \code{lambda}.
#' Fits Gaussian (linear), binomial Logit (logistic), Poisson, multinomial Logit regression models, and
#' Cox proportional hazard model with various non-convex penalties.
#' @param y.vec (numeric vector) response vector.
#' Must be 0,1 for \code{binomial} and 1,2,..., for \code{multinomial}.
#' @param x.mat (numeric matrix) design matrix without intercept.
#' The censoring indicator must be included at the last column of the design matrix for \code{cox}.
#' @param family (character) regression model. Supported models are
#' \code{gaussian},
#' \code{binomial},
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
#' @param add.max (numeric) the maximum number of variables added in CCCP iterations (largely internal use. See referecnes).
#' @param niter.max (numeric) maximum number of iterations in CCCP.
#' @param qiter.max (numeric) maximum number of quadratic approximations in each CCCP iteration.
#' @param aiter.max (numeric) maximum number of iterations in CD algorithm.
#' @param b.eps (numeric) convergence threshold for coefficients vector.
#' @param k.eps (numeric) convergence threshold for KKT conditions.
#' @param c.eps (numeric) convergence threshold for KKT conditions (largely internal use).
#' @param cut (logical) convergence threshold for KKT conditions  (largely internal use).
#' @param local (logical) whether to use local initial estimator for path construction. It may take a long time.
#' @param local.initial (nemeric vector) initial estiamtor for \code{local=TRUE}.
#' @details
#' The sequence of models indexed by \code{lambda} is fit
#' by using concave convex procedure (CCCP) and coordinate descent (CD) algorithm (see references).
#' The objective function is \deqn{ (sum of squared residuals)/2n + [alpha*penalty + (1-alpha)*ridge] } for \code{gaussian}
#' and \deqn{ (log-likelihood)/n - [alpha*penalty + (1-alpha)*ridge] } for the others,
#' assuming the canonical link.
#' The algorithm applies the warm start strategy (see references) and tries projections
#' after \code{proj.min} iterations in CD algorithm, which makes the algorithm fast and stable.
#' \code{x.standardize} makes each column of \code{x.mat} to have the same Euclidean length
#' but the coefficients will be rescaled into the original.
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
#'   \item{beta}{fitted coefficients. Use \code{coef.ncpen} for \code{multinomial} since the coefficients are representd as vectors.}
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
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="binomial",penalty="classo")
#' ### poison regression with mlog penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="poisson")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="poisson",penalty="mlog")
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="multinomial",penalty="sridge")
#' coef.ncpen(fit)
#' ### cox regression with mbridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,r=0.2,cf.min=0.5,cf.max=1,corr=0.5,family="cox")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="cox",penalty="mbridge")
#' @export
ncpen = function(y.vec,x.mat,
                 family=c("gaussian","binomial","multinomial","cox","poisson"),
                 penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                 x.standardize=TRUE,intercept=TRUE,
                 lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                 df.max=50,cf.max=100,proj.min=10,add.max=10,niter.max=30,qiter.max=10,aiter.max=100,
                 b.eps=1e-7,k.eps=1e-4,c.eps=1e-6,cut=TRUE,local=FALSE,local.initial=NULL){
  family = match.arg(family); penalty = match.arg(penalty)
  tun = control.ncpen(y.vec,x.mat,family,penalty,x.standardize,intercept,
                  lambda,n.lambda,r.lambda,w.lambda,gamma,tau,alpha,aiter.max,b.eps)
  if(local==TRUE){
    if(is.null(local.initial)){ stop(" supply 'local.initial' \n") }
    warning("       'local==TRUE' option may take a long time \n")
  } else { local.initial = rep(0,length(tun$w.lambda)) }
  fit = native_cpp_ncpen_fun_(tun$y.vec,tun$x.mat,
                              tun$w.lambda,tun$lambda,tun$gamma,tun$tau,tun$alpha,
                              df.max,niter.max,qiter.max,aiter.max,
                              b.eps,k.eps,proj.min,cut,c.eps,add.max,family,penalty,local,local.initial,cf.max)
  if(x.standardize==TRUE){ fit$beta = fit$beta/tun$std  }
  ret = list(y.vec=tun$y.vec,x.mat=tun$x.mat,
             family=family,penalty=penalty,
             x.standardize=x.standardize,intercept=tun$intercept,std=tun$std,
             lambda=drop(fit$lambda),w.lambda=tun$w.lambda,gamma=tun$gamma,tau=tun$tau,alpha=tun$alpha,
             local.initial=local.initial,
             beta=fit$beta,df=drop(fit$df))

  class(ret) = "ncpen";
  return(ret);
}
######################################################### cv.ncpen
#' @title
#' cv.ncpen: cross validation for \code{ncpen}
#' @description
#' performs k-fold cross-validation (CV) for nonconvex penalized regression models
#' over a sequence of the regularization parameter \code{lambda}.
#' @param y.vec (numeric vector) response vector.
#' Must be 0,1 for \code{binomial} and 1,2,..., for \code{multinomial}.
#' @param x.mat (numeric matrix) design matrix without intercept.
#' The censoring indicator must be included at the last column of the design matrix for \code{cox}.
#' @param n.fold (numeric) number of folds for CV.
#' @param fold.id (numeric vector) fold ids from 1 to k that indicate fold configuration.
#' @param family (character) regression model. Supported models are
#' \code{gaussian},
#' \code{binomial},
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
#' @param add.max (numeric) the maximum number of variables added in CCCP iterations (largely internal use. See referecnes).
#' @param niter.max (numeric) maximum number of iterations in CCCP.
#' @param qiter.max (numeric) maximum number of quadratic approximations in each CCCP iteration.
#' @param aiter.max (numeric) maximum number of iterations in CD algorithm.
#' @param b.eps (numeric) convergence threshold for coefficients vector.
#' @param k.eps (numeric) convergence threshold for KKT conditions.
#' @param c.eps (numeric) convergence threshold for KKT conditions (largely internal use).
#' @param cut (logical) convergence threshold for KKT conditions  (largely internal use).
#' @param local (logical) whether to use local initial estimator for path construction. It may take a long time.
#' @param local.initial (nemeric vector) initial estiamtor for \code{local=TRUE}.
#' @details
#' Two kinds of CV errors are returned: root mean squard error and negative log likelihood.
#' The resutls depends on the random partition made interally.
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
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10)
#' coef.cv.ncpen(fit)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="binomial",penalty="classo")
#' coef.cv.ncpen(fit)
#' ### poison regression with mlog penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="poisson")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="poisson",penalty="mlog")
#' coef.cv.ncpen(fit)
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="multinomial",penalty="sridge")
#' coef.cv.ncpen(fit)
#' ### cox regression with mcp penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,r=0.2,cf.min=0.5,cf.max=1,corr=0.5,family="cox")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="cox",penalty="scad")
#' coef.cv.ncpen(fit)
#' @export
cv.ncpen = function(y.vec,x.mat,n.fold=10,fold.id=NULL,
                    family=c("gaussian","binomial","multinomial","cox","poisson"),
                    penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                    x.standardize=TRUE,intercept=TRUE,
                    lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                    df.max=50,cf.max=100,proj.min=10,add.max=10,niter.max=30,qiter.max=10,aiter.max=100,
                    b.eps=1e-6,k.eps=1e-4,c.eps=1e-6,cut=TRUE,local=FALSE,local.initial=NULL){
  family = match.arg(family); penalty = match.arg(penalty)
  fit = ncpen(y.vec,x.mat,family,penalty,
              x.standardize,intercept,
              lambda,n.lambda,r.lambda,w.lambda,gamma,tau,alpha,
              df.max,cf.max,proj.min,add.max,niter.max,qiter.max,aiter.max,
              b.eps,k.eps,c.eps,cut,local,local.initial)
  if(is.null(fold.id)){
    if(family=="cox"){ c.vec = x.mat[,dim(x.mat)[2]] } else { c.vec = y.vec }
    idx = fold.cv.ncpen(c.vec,n.fold,family)$idx
  } else { idx = fold.id }
  if(min(table(idx))<4) stop("at least 3 samples in each fold")
  rmse = like = list()
  for(s in 1:n.fold){
    cat("cv fold number:",s,"\n")
    yset = idx==s; xset = yset; if(family=="multinomial"){ oxset = xset; xset = rep(xset,max(fit$y.vec)-1) }
    fold.fit = native_cpp_ncpen_fun_(fit$y.vec[!yset],fit$x.mat[!xset,],
                                     fit$w.lambda,fit$lambda,fit$gamma,fit$tau,fit$alpha,
                                     df.max,niter.max,qiter.max,aiter.max,
                                     b.eps,k.eps,proj.min,cut,c.eps,add.max,family,penalty,
                                     local,fit$local.initial,cf.max)
    if(x.standardize==TRUE){ fold.fit$beta = fold.fit$beta/fit$std }
    fold.fit$intercept = fit$intercept; fold.fit$family = family; fold.fit$y.vec = fit$y.vec[!yset]
    if(family!="multinomial"){
      if(fit$int==TRUE){
        rmse[[s]] = predict.ncpen(fold.fit,type="rmse",new.y.vec=fit$y.vec[yset],new.x.mat=fit$x.mat[xset,-1])
        like[[s]] = predict.ncpen(fold.fit,type="like",new.y.vec=fit$y.vec[yset],new.x.mat=fit$x.mat[xset,-1])
      } else {
        rmse[[s]] = predict.ncpen(fold.fit,type="rmse",new.y.vec=fit$y.vec[yset],new.x.mat=fit$x.mat[xset,])
        like[[s]] = predict.ncpen(fold.fit,type="like",new.y.vec=fit$y.vec[yset],new.x.mat=fit$x.mat[xset,])
      }
    } else {
      if(fit$int==TRUE){
        rmse[[s]] = predict.ncpen(fold.fit,type="rmse",new.y.vec=fit$y.vec[yset],
                                  new.x.mat=fit$x.mat[1:length(oxset),][oxset,2:(dim(x.mat)[2]+1)])
        like[[s]] = predict.ncpen(fold.fit,type="like",new.y.vec=fit$y.vec[yset],
                                  new.x.mat=fit$x.mat[1:length(oxset),][oxset,2:(dim(x.mat)[2]+1)])
      } else {
        rmse[[s]] = predict.ncpen(fold.fit,type="rmse",new.y.vec=fit$y.vec[yset],
                                  new.x.mat=fit$x.mat[1:length(oxset),][oxset,1:dim(x.mat)[2]])
        like[[s]] = predict.ncpen(fold.fit,type="like",new.y.vec=fit$y.vec[yset],
                                  new.x.mat=fit$x.mat[1:length(oxset),][oxset,1:dim(x.mat)[2]])
      }
    }
  }
  nlam = min(unlist(lapply(rmse,length)))
  rmse = rowSums(matrix(unlist(lapply(rmse,function(x) x[1:nlam])),ncol=n.fold))
  nlam = min(unlist(lapply(like,length)))
  like = rowSums(matrix(unlist(lapply(like,function(x) x[1:nlam])),ncol=n.fold))
  ret = list(ncpen.fit=fit,
             fold.index = idx,
             rmse = rmse,
             like = like,
             lambda = fit$lambda[1:nlam]
  )
  class(ret) = "cv.ncpen"
  return(ret)
}
######################################################### control.ncpen
#' @title
#' control.ncpen: do preliminary works for \code{ncpen}.
#' @description
#' The function returns controlled samples and tuning parameters for \code{ncpen} by eliminating unnecessary errors.
#' @param y.vec (numeric vector) response vector.
#' Must be 0,1 for \code{binomial} and 1,2,..., for \code{multinomial}.
#' @param x.mat (numeric matrix) design matrix without intercept.
#' The censoring indicator must be included at the last column of the design matrix for \code{cox}.
#' @param family (character) regression model. Supported models are
#' \code{gaussian},
#' \code{binomial},
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
#' @param aiter.max (numeric) maximum number of iterations in CD algorithm.
#' @param b.eps (numeric) convergence threshold for coefficients vector.

#' @details
#' The function is used internal purpose but useful when users want to extract proper tuning parameters for \code{ncpen}.
#' Do not supply the samples from \code{control.ncpen} into \code{ncpen} or \code{cv.ncpen} directly to avoid unexpected errors.
#' @return An object with S3 class \code{ncpen}.
#'   \item{y.vec}{response vector.}
#'   \item{x.mat}{design matrix adjusted to supplied options such as family and intercept.}
#'   \item{family}{regression model.}
#'   \item{penalty}{penalty.}
#'   \item{x.standardize}{whether to standardize \code{x.mat}.}
#'   \item{intercept}{whether to include the intercept.}
#'   \item{std}{scale factor for \code{x.standardize}.}
#'   \item{lambda}{lambda values for the analysis.}
#'   \item{n.lambda}{the number of \code{lambda} values.}
#'   \item{r.lambda}{ratio of the smallest \code{lambda} value to largest.}
#'   \item{w.lambda}{penalty weights for each coefficient.}
#'   \item{gamma}{additional tuning parameter for controlling shrinkage effect of \code{classo} and \code{sridge} (see references).}
#'   \item{tau}{concavity parameter of the penalties (see references).}
#'   \item{alpha}{ridge effect (amount of ridge penalty). see details.}
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
#' \code{\link{ncpen}}, \code{\link{cv.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' tun = control.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,tau=1)
#' tun$tau
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' tun = control.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,
#'                     family="multinomial",penalty="sridge",gamma=10)
#' ### cox regression with mcp penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,r=0.2,cf.min=0.5,cf.max=1,corr=0.5,family="cox")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' tun = control.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="cox",penalty="scad")
#' @export
control.ncpen = function(y.vec,x.mat,
                         family=c("gaussian","binomial","multinomial","cox","poisson"),
                         penalty=c("scad","mcp","tlp","lasso","classo","ridge","sridge","mbridge","mlog"),
                         x.standardize=TRUE,intercept=TRUE,
                         lambda=NULL,n.lambda=NULL,r.lambda=NULL,w.lambda=NULL,gamma=NULL,tau=NULL,alpha=NULL,
                         aiter.max=1e+2,b.eps = 1e-7){
  family = match.arg(family); penalty = match.arg(penalty);
  if(is.null(alpha)){ alpha = switch(penalty,scad=1,mcp=1,tlp=1,lasso=1,classo=1,ridge=0,sridge=1,mbridge=1,mlog=1)
  } else {
    if(penalty=="ridge"){ if(alpha!=0){ warning("alpha is set to 0 \n") }; alpha = 0;
    } else if(penalty=="sridge"){ if(alpha!=1){ warning("alpha is set to 1 \n") }; alpha = 1;
    } else { if(alpha<0){ warning("alpha is set to 0 \n"); alpha = 0 }; if(alpha>1){ warning("alpha is set to 1 \n"); alpha = 1 } }
  }
  if(is.null(tau)){ tau = switch(penalty,scad=3.7,mcp=2.1,tlp=0.001,lasso=2,classo=2.1,ridge=2,sridge=2.1,mbridge=0.001,mlog=0.001)
  } else {
    if(penalty=="scad"){ if(tau<=2){ tau = 2.001; warning("tau is set to ",tau,"\n") }
    } else if(penalty=="mcp"){ if(tau<=1){ tau = 1.001; warning("tau is set to ",tau,"\n") }
    } else if((penalty=="classo")|(penalty=="sridge")){ if(tau<=1){ tau = 1.001; warning("tau is set to ",tau,"\n") }
    } else { if(tau<=0){ tau = 0.001; warning("tau is set to ",tau,"\n") } }
  }
  n = dim(x.mat)[1]; p = dim(x.mat)[2]; if(family=="cox"){ p = p-1 }
  if(is.null(w.lambda)){ w.lambda = rep(1,p)
  } else {
    if(length(w.lambda)!=p){ stop("the number of elements in w.lambda should be the number of input variables")
    } else if(sum(w.lambda==0)>=n){ stop("the number of zero elements in w.lambda should be less than the sample size")
    } else if(min(w.lambda)<0){ stop("elements in w.lambda should be non-negative")
    } else { w.lambda = (p-sum(w.lambda==0))*w.lambda/sum(w.lambda) }
  }
  if(x.standardize==FALSE){ std = rep(1,p)
  } else {
    if(family=="cox"){ std = sqrt(colSums(x.mat[,-(p+1)]^2)/n); std[std==0] = 1; x.mat[,-(p+1)] = sweep(x.mat[,-(p+1)],2,std,"/")
    } else { std = sqrt(colSums(x.mat^2)/n); std[std==0] = 1; x.mat = sweep(x.mat,2,std,"/") }
  }
  if(family=="cox") intercept = FALSE;
  if(intercept==TRUE){ x.mat = cbind(1,x.mat); std = c(1,std); w.lambda = c(0,w.lambda); p = p+1; }
  if(family=="multinomial"){
    k = max(y.vec)
    if(length(unique(y.vec))!=k){ stop("label must be denoted by 1,2,...: ") }
    x.mat = kronecker(diag(k-1),x.mat); w.lambda = rep(w.lambda,k-1); p = p*(k-1); std = rep(std,k-1);
  }
  if(is.null(r.lambda)){
    if((family=="cox")|(family=="multinomial")){
      r.lambda = ifelse(n<p,0.1,0.01);
    } else {
      r.lambda = ifelse(n<p,0.01,0.001);
    }
  } else {
    if((r.lambda<=0)|(r.lambda>=1)) stop("r.lambda should be between 0 and 1");
  }
  if(is.null(n.lambda)){ n.lambda = 100;
  } else { if((n.lambda<1e+1)|(n.lambda>1e+3)) stop("n.lambda should be between 10 and 1000") }
  if(sum(w.lambda==0)==0){ b.vec = rep(0,p);
  } else { b.vec = rep(0,p);
  if(family!="cox"){ a.set = c(1:p)[w.lambda==0]; } else { a.set = c(c(1:p)[w.lambda==0],p+1); }
  ax.mat = x.mat[,a.set,drop=F]; b.vec[a.set] = drop(native_cpp_nr_fun_(family,y.vec,ax.mat,aiter.max,b.eps));
  }
  g.vec = abs(native_cpp_obj_grad_fun_(family,y.vec,x.mat,b.vec))
  if(alpha!=0){
    lam.max = max(g.vec[w.lambda!=0]/w.lambda[w.lambda!=0])/alpha; lam.max = lam.max+lam.max/10
    lam.vec = exp(seq(log(lam.max),log(lam.max*r.lambda),length.out=n.lambda))
    if(is.null(lambda)){ lambda = lam.vec
    } else {
      if(min(lambda)<=0){ stop("elements in lambda should be positive");
      } else { lambda = sort(lambda,decreasing=TRUE);
      if(max(lambda)<lam.max){ lambda = c(lam.vec[lam.vec>lambda[1]],lambda); warning("lambda is extended up to ",lam.max,"\n") }
      }
    }
  }
  if(alpha==0) lambda = exp(seq(log(1e-2),log(1e-5),length.out=n.lambda))
  if(is.null(gamma)){
    gamma = switch(penalty,scad=0,mcp=0,tlp=0,lasso=0,classo=min(lambda)/2,ridge=0,sridge=min(lambda)/2,mbridge=0,mlog=0)
  } else {
    if((penalty!="classo")&(penalty!="sridge")){ gamma = 0
    } else {
      if(gamma<0){ gamma = min(lambda)/2; warning("gamma is negative and set to ",gamma,"\n") }
      if(gamma>=max(lambda)){ gamma = max(lambda)/2; warning("gamma is larger than lambda and set to ",gamma,"\n") }
      if(gamma>min(lambda)){ lambda = lambda[lambda>=gamma]; warning("lambda is set larger than gamma=",gamma,"\n") }
    }
  }
  ret = list(y.vec=y.vec,x.mat=x.mat,
             family=family,penalty=penalty,
             x.standardize=x.standardize,intercept=intercept,std=std,
             lambda=lambda,n.lambda=n.lambda,r.lambda=r.lambda,w.lambda=w.lambda,gamma=gamma,tau=tau,alpha=alpha)
  class(ret) = "control.ncpen";
  return(ret)
}
######################################################### coef.ncpen
#' @title
#' coef.ncpen: extract the coefficients from an \code{ncpen} object
#' @description
#' The function returns the coefficients matrix for all lambda values.
#' @param object (ncpen object) fitted \code{ncpen} object.
#' @param ... other S3 parameters. Not used.
#' @return
#'  \item{beta}{The coefficients matrix or list for \code{multinomial}.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat)
#' coef(fit)
#' ### multinomial regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="multinomial",penalty="classo")
#' coef(fit)
#' @export coef.ncpen
#' @export
coef.ncpen = function(object, ...){
     fit = object;
  p = dim(fit$beta)[1]
  k = max(fit$y.vec)
  p = ifelse(fit$fam=="multinomial",p/(k-1),p)
  p = ifelse(fit$int==T,p-1,p)
  lam.name = paste("lambda",1:length(fit$lam),sep="")
  var.name = paste("x",1:p,sep="")
  if(fit$int==TRUE){ var.name = c("intercept",var.name) }
  beta = fit$beta
  if(fit$fam=="multinomial"){
    beta = list();
    for(i in 1:dim(fit$beta)[2]){
      beta[[i]] = cbind(matrix(fit$beta[,i],ncol=k-1),0)
      rownames(beta[[i]]) = var.name
      colnames(beta[[i]]) = paste("class",1:k,sep="")
    }
    names(beta) = lam.name
  } else {
    colnames(beta) = lam.name; rownames(beta) = var.name
  }
  return(beta)
}
######################################################### plot.ncpen
#' @title
#' plot.ncpen: plots coefficients from an \code{ncpen} object.
#' @description
#' Produces a plot of the coefficients paths for a fitted \code{ncpen} object. Class-wise paths can be drawn for \code{multinomial}.
#' @param x (ncpen object) Fitted \code{ncpen} object.
#' @param log.scale (logical) whether to use log scale of lambda for horizontal axis.
#' @param mult.type (character) additional option for \code{multinomial} whether to draw the coefficients class-wise or not.
#' Default is \code{mat} that uses class-wise coefficients.
#' @param ... other graphical parameters to \code{\link{plot}}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat)
#' plot(fit)
#' ### multinomial regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="multinomial",penalty="classo")
#' plot(fit)
#' plot(fit,mult.type="vec",log.scale=TRUE)
#' @export plot.ncpen
#' @export
plot.ncpen = function(x,log.scale=FALSE,mult.type=c("mat","vec"),...){
  fit = x;
  beta = fit$beta
  lambda = fit$lam
  if(log.scale==T) lambda = log(lambda)
  if(fit$fam!="multinomial"){
    if(fit$int==T){ beta = beta[-1,] }
    plot(lambda,beta[1,],type="n",ylim=c(min(beta),max(beta)),main="trace of coefficients",ylab="",...)
    for(j in 1:dim(beta)[1]){ lines(lambda,beta[j,],col=j,...) } #col=iter,lty=iter) }
  } else {
    k = max(fit$y.vec)
    p = dim(beta)[1]/(k-1)
    if(match.arg(mult.type)=="mat"){#"mat
      par(mfrow=c(1,k-1))
      for(i in 1:(k-1)){
        mat = beta[(1+(i-1)*p):(i*p),]
        if(fit$int==T){ mat = mat[-1,] }
        plot(lambda,mat[1,],type="n",ylim=c(min(mat),max(mat)),main="trace of coefficients",ylab="",...)
        for(j in 1:dim(mat)[1]){ lines(lambda,mat[j,],col=j,...) } #col=iter,lty=iter) }
      }
    } else {#"vec
      if(fit$int==T){ beta = beta[-(1+(0:(k-2))*p),] }
      plot(lambda,beta[1,],type="n",ylim=c(min(beta),max(beta)),main="trace of coefficients",ylab="",...)
      for(j in 1:dim(beta)[1]){ lines(lambda,beta[j,],col=j,...) } #col=iter,lty=iter) }
    }
  }
}
######################################################### gic.ncpen
#' @title
#' gic.ncpen: compute the generalized information criterion (GIC) for the selection of lambda
#' @description
#' The function provides the selection of the regularization parameter lambda based
#' on the GIC including AIC and BIC.
#' @param fit (ncpen object) fitted \code{ncpen} object.
#' @param weight (numeric) the weight factor for various information criteria.
#' Default is BIC if \code{n>p} and GIC if \code{n<p} (see details).
#' @param verbose (logical) whether to plot the GIC curve.
#' @param ... other graphical parameters to \code{\link{plot}}.
#' @details
#' User can supply various \code{weight} values (see references). For example,
#' \code{weight=2},
#' \code{weight=log(n)},
#' \code{weight=log(log(p))log(n)},
#' \code{weight=log(log(n))log(p)},
#' corresponds to AIC, BIC (fixed dimensional model), modified BIC (diverging dimensional model) and GIC (high dimensional model).
#' @return The coefficients \code{\link{matrix}}.
#'   \item{gic}{the GIC values.}
#'   \item{lambda}{the sequence of lambda values used to calculate GIC.}
#'   \item{opt.beta}{the optimal coefficients selected by GIC.}
#'   \item{opt.lambda}{the optimal lambda value.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Wang, H., Li, R. and Tsai, C.L. (2007). Tuning parameter selectors for the smoothly clipped absolute deviation method.
#' \emph{Biometrika}, 94(3), 553-568.
#' Wang, H., Li, B. and Leng, C. (2009). Shrinkage tuning parameter selection with a diverging number of parameters.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(3), 671-683.
#' Kim, Y., Kwon, S. and Choi, H. (2012). Consistent Model Selection Criteria on High Dimensions.
#' \emph{Journal of Machine Learning Research}, 13, 1037-1057.
#' Fan, Y. and Tang, C.Y. (2013). Tuning parameter selection in high dimensional penalized likelihood.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 75(3), 531-552.
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat)
#' gic.ncpen(fit,pch="*",type="b")
#' ### multinomial regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec,x.mat=x.mat,family="multinomial",penalty="classo")
#' gic.ncpen(fit,pch="*",type="b")
#' @export gic.ncpen
#' @export
gic.ncpen = function(fit,weight=NULL,verbose=TRUE,...){
  n = length(fit$y.vec)
  p = dim(fit$beta)[1]
  k = max(fit$y.vec)
  if(is.null(weight)){ weight = ifelse(n>p,log(n),log(log(n))*log(p)) }
  dev = 2*n*apply(fit$beta,2,FUN="native_cpp_obj_fun_",name=fit$fam,y_vec=fit$y.vec,x_mat=fit$x.mat)
  gic = dev+weight*drop(fit$df)
  opt = which.min(gic)
  if(verbose==TRUE){
    plot(fit$lam,gic,xlab="lambda",ylab="",main="information criterion",...)
    abline(v=fit$lam[opt])
  }
  opt.beta=fit$beta[,opt]
  if(fit$fam=="multinomial"){
    opt.beta = cbind(matrix(fit$beta[,opt],ncol=k-1),0); colnames(opt.beta) = paste("class",1:k,sep="")
  }
  return(list(gic=gic,lambda=fit$lambda,opt.lambda=fit$lambda[opt],opt.beta=opt.beta))
}
######################################################### predict.ncpen
#' @title
#' predict.ncpen: make predictions from an \code{ncpen} object
#' @description
#' The function provides various types of predictions from a fitted \code{ncpen} object:
#' response, regression, probability, root mean squared error (RMSE), negative log-likelihood (LIKE).
#' @param object (ncpen object) fitted \code{ncpen} object.
#' @param new.y.vec (numeric vector). vector of new response at which predictions are to be made.
#' @param new.x.mat (numeric matrix). matrix of new design at which predictions are to be made.
#' @param type (character) type of prediction.
#' \code{y} returns new responses from \code{new.x.mat}.
#' \code{reg} returns new linear predictors from \code{new.x.mat}.
#' \code{prob} returns new class probabilities from \code{new.x.mat} for \code{binomial} and \code{multinomial}.
#' \code{rmse} returns RMSE from \code{new.y.vec} and \code{new.x.mat}.
#' \code{prob} returns LIKE from \code{new.y.vec} and \code{new.x.mat}.
#' @param prob.cut (numeric) threshold value of probability for \code{binomial}.
#' @param ... other S3 parameters. Not used.
#' @return prediction values depending on \code{type} for all lambda values.
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec[1:190],x.mat=x.mat[1:190,])
#' predict(fit,"y",new.x.mat=x.mat[190:200,])
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec[1:190],x.mat=x.mat[1:190,],family="binomial",penalty="classo")
#' predict(fit,"y",new.x.mat=x.mat[190:200,])
#' predict(fit,"y",new.x.mat=x.mat[190:200,],prob.cut=0.3)
#' predict(fit,"reg",new.x.mat=x.mat[190:200,])
#' predict(fit,"prob",new.x.mat=x.mat[190:200,])
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = ncpen(y.vec=y.vec[1:190],x.mat=x.mat[1:190,],family="multinomial",penalty="classo")
#' predict(fit,"y",new.x.mat=x.mat[190:200,])
#' predict(fit,"reg",new.x.mat=x.mat[190:200,])
#' predict(fit,"prob",new.x.mat=x.mat[190:200,])
#' @export predict.ncpen
#' @export
predict.ncpen = function(object,type=c("y","reg","prob","rmse","like"),new.y.vec=NULL,new.x.mat=NULL,prob.cut=0.5, ...){
     fit = object;

  if(is.null(new.x.mat)){ stop("'new.x.mat' should be supplied for prediction") }
  if(match.arg(type)%in%c("y","reg","prob")){ new.y.vec = NULL
  } else { if(is.null(new.y.vec)) stop("'new.y.vec' should be supplied for prediction") }
  l.name = paste("lam",1:length(fit$lam),sep="")
  if(fit$int==TRUE){ new.x.mat = cbind(1,new.x.mat) }
  if(fit$fam=="multinomial"){
    k = max(fit$y.vec); new.x.mat = kronecker(diag(k-1),new.x.mat); c.name = paste("class",1:k,sep="")
  }
  if(fit$fam!="cox"){
    new.reg = new.x.mat%*%fit$beta
  } else { new.reg = new.x.mat[,-dim(new.x.mat)[2]]%*%fit$beta }

  if(fit$fam=="gaussian"){ new.y = new.reg; new.prob = NULL
  } else if(fit$fam=="binomial"){ new.prob = exp(new.reg)/(1+exp(new.reg)); new.y = 1*(new.prob>prob.cut)
  } else if(fit$fam=="multinomial"){
    reg.list = list()
    for(i in 1:dim(new.reg)[2]){ reg.list[[i]] = cbind(matrix(new.reg[,i],ncol=k-1),0); colnames(reg.list[[i]]) = c.name }
    new.reg = reg.list; new.prob = lapply(reg.list,function(x) exp(x)/rowSums(exp(x)))
    new.y = matrix(unlist(lapply(new.prob,function(x) apply(x,1,which.max))),ncol=length(fit$lam))
  } else if(fit$fam=="poisson"){ new.y = exp(new.reg); new.prob = NULL
  } else { new.y = exp(new.reg); new.prob = NULL }

  if(match.arg(type)=="reg"){ return(reg=new.reg)
  } else if(match.arg(type)=="prob"){ return(prob=new.prob)
  } else if(match.arg(type)=="y"){ return(y=new.y)
  } else if(match.arg(type)=="rmse"){
    new.rmse = sqrt(colSums((new.y.vec-new.y)^2)); return(rmse=new.rmse)
  } else {
    new.like = apply(fit$beta,2,FUN="native_cpp_obj_fun_",name=fit$fam,y_vec=new.y.vec,x_mat=new.x.mat)
    return(rmse=new.like)
  }
}
######################################################### fold.cv.ncpen
#' @title
#' fold.cv.ncpen: extracts fold ids for \code{cv.ncpen}.
#' @description
#' The function returns fold configuration of the samples for CV.
#' @param c.vec (numeric vector) vector for construction of CV ids:
#' censoring indicaor for \code{cox} and response vector for the others.
#' @param n.fold (numeric) number of folds for CV.
#' @param family (character) regression model. Supported models are
#' \code{gaussian},
#' \code{binomial},
#' \code{poisson},
#' \code{multinomial},
#' and \code{cox}.
#' Default is \code{gaussian}.
#' @return fold ids of the samples.
#'  \item{idx}{fold ids.}
#'  \item{n.fold}{the number of folds.}
#'  \item{family}{the model.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{cv.ncpen}}, \code{\link{plot.cv.ncpen}} , \code{\link{gic.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fold.id = fold.cv.ncpen(c.vec=y.vec,n.fold=10)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fold.id = fold.cv.ncpen(c.vec=y.vec,n.fold=10,family="binomial")
#' ### poison regression with mlog penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="poisson")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fold.id = fold.cv.ncpen(c.vec=y.vec,n.fold=10,family="poisson")
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fold.id = fold.cv.ncpen(c.vec=y.vec,n.fold=10,family="multinomial")
#' ### cox regression with mcp penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,r=0.2,cf.min=0.5,cf.max=1,corr=0.5,family="cox")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fold.id = fold.cv.ncpen(c.vec=x.mat[,21],n.fold=10,family="cox")
#' @export
fold.cv.ncpen = function(c.vec,n.fold=10,family=c("gaussian","binomial","multinomial","cox","poisson")){
  family = match.arg(family)
  n = length(c.vec); idx = c.vec*0
  if(family=="gaussian"){ c.vec = c.vec-median(c.vec)+0.5 }
  if(family!="multinomial"){
    set0 = c(1:n)[c.vec<0.5]; idx[set0] = sample(rep(1:n.fold,length=length(set0)))
    set1 = c(1:n)[-set0]; idx[set1] = sample(rep(1:n.fold,length=length(set1)))
  } else { u.vec = unique(c.vec); k = length(u.vec);
  for(i in 1:k){ set = c(1:n)[c.vec==u.vec[i]]; idx[set] = sample(rep(1:n.fold,length=length(set))) }
  }
  ret = list(family=family,n.fold=n.fold,idx=idx)
  return(ret)
}
######################################################### coef.cv.ncpen
#' @title
#' coef.cv.ncpen: extracts the optimal coefficients from \code{cv.ncpen}.
#' @description
#' The function returns the optimal vector of coefficients.
#' @param object (cv.ncpen object) fitted \code{cv.ncpen} object.
#' @param type (character) a cross-validated error type which is either \code{rmse} or \code{like}.
#' @param ... other S3 parameters. Not used.
#' Each error type is defined in \code{\link{cv.ncpen}}.
#' @return the optimal coefficients vector selected by cross-validation.
#'  \item{type}{error type.}
#'  \item{lambda}{the optimal lambda selected by  CV.}
#'  \item{beta}{the optimal coefficients selected by CV.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{cv.ncpen}}, \code{\link{plot.cv.ncpen}} , \code{\link{gic.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10)
#' coef(fit)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="binomial",penalty="classo")
#' coef(fit)
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=10,family="multinomial",penalty="sridge")
#' coef(fit)
#' @export coef.cv.ncpen
#' @export
coef.cv.ncpen = function(object,type=c("rmse","like"), ...){
     cv.fit = object;

  type = match.arg(type)
  if(type=="rmse"){ opt = which.min(cv.fit$rmse) }
  if(type=="like"){ opt = which.min(cv.fit$like) }
  beta=cv.fit$ncpen.fit$beta[,opt]
  if(cv.fit$ncpen.fit$fam=="multinomial"){ k = max(cv.fit$ncpen.fit$y.vec); beta = cbind(matrix(beta,ncol=k-1),0) }
  return(list(type=match.arg(type),lambda=cv.fit$ncpen.fit$lambda[opt],beta=beta))
}
######################################################### plot.cv.ncpen
#' @title
#' plot.cv.ncpen: plot cross-validation error curve.
#' @description
#' The function Produces a plot of the cross-validated errors from \code{cv.ncpen} object.
#' @param x fitted \code{cv.ncpen} object.
#' @param type (character) a cross-validated error type which is either \code{rmse} or \code{like}.
#' @param log.scale (logical) whether to use log scale of lambda for horizontal axis.
#' @param ... other graphical parameters to \code{\link{plot}}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Lee, S., Kwon, S. and Kim, Y. (2016). A modified local quadratic approximation algorithm for penalized optimization problems.
#' \emph{Computational Statistics and Data Analysis}, 94, 275-286.
#' @seealso
#' \code{\link{cv.ncpen}}
#' @examples
#' ### linear regression with scad penalty
#' par(mfrow=c(1,2))
#' sam =  sam.gen.ncpen(n=500,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=50)
#' plot(fit)
#' plot(fit,log.scale=F)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=50,family="binomial",penalty="scad")
#' plot(fit,type="rmse")
#' plot(fit,type="like")
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=10,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' fit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,n.lambda=50,family="multinomial",penalty="sridge")
#' plot(fit,type="like")
#' @export plot.cv.ncpen
#' @export
plot.cv.ncpen = function(x,type=c("rmse","like"),log.scale=FALSE,...){
     cv.fit = x;

     type = match.arg(type)
     if(type=="rmse"){
          opt = which.min(cv.fit$rmse)
          plot(cv.fit$lambda,cv.fit$rmse,main="root mean square error",type="b",pch="*")
     }
     if(type=="like"){
          opt = which.min(cv.fit$like)
          plot(cv.fit$lambda,cv.fit$like,main="negative log-likelihood",type="b",pch="*")
     }
     abline(v=cv.fit$lambda[opt])
}
######################################################### sam.gen.ncpen
#' @title
#' sam.gen.ncpen: generate a simulated dataset.
#' @description
#' Generate a synthetic dataset based on the correlation structure from generalized linear models.
#' @param n (numeric) the number of samples.
#' @param p (numeric) the number of variables.
#' @param q (numeric) the number of nonzero coefficients.
#' @param k (numeric) the number of classes for \code{multinomial}.
#' @param r (numeric) the ratio of censoring for \code{cox}.
#' @param cf.min (numeric) value of the minimum coefficient.
#' @param cf.max (numeric) value of the maximum coefficient.
#' @param corr (numeric) strength of correlations in the correlation structure.
#' @param family (character) model type.
#' @param seed (numeric) seed number for random generation. Default does not use seed.
#' @details
#' A design matrix for regression models is generated from the multivariate normal distribution with a correlation structure.
#' Then the response variables are computed with a specific model based on the true coefficients (see references).
#' Note the censoring indicator locates at the last column of \code{x.mat} for \code{cox}.
#' @return An object with list class containing
#'   \item{x.mat}{design matrix.}
#'   \item{y.vec}{responses.}
#'   \item{b.vec}{true coefficients.}
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#' @references
#' Kwon, S., Lee, S. and Kim, Y. (2016). Moderately clipped LASSO.
#' \emph{Computational Statistics and Data Analysis}, 92C, 53-67.
#' Kwon, S. and Kim, Y. (2012). Large sample properties of the SCAD-penalized maximum likelihood estimation on high dimensions.
#' \emph{Statistica Sinica}, 629-653.
#' @seealso
#' \code{\link{ncpen}}
#' @examples
#' ### linear regression
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5)
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' head(x.mat); head(y.vec)
#' ### logistic regression with classo penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="binomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' head(x.mat); head(y.vec)
#' ### poison regression with mlog penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,cf.min=0.5,cf.max=1,corr=0.5,family="poisson")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' head(x.mat); head(y.vec)
#' ### multinomial regression with sridge penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,k=3,cf.min=0.5,cf.max=1,corr=0.5,family="multinomial")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' head(x.mat); head(y.vec)
#' ### cox regression with mcp penalty
#' sam =  sam.gen.ncpen(n=200,p=20,q=5,r=0.2,cf.min=0.5,cf.max=1,corr=0.5,family="cox")
#' x.mat = sam$x.mat; y.vec = sam$y.vec
#' head(x.mat); head(y.vec)
#' @export
sam.gen.ncpen = function(n=100,p=50,q=10,k=3,r=0.3,cf.min=0.5,cf.max=1,corr=0.5,seed=NULL,
                         family=c("gaussian","binomial","multinomial","cox","poisson")){
  if(!is.null(seed)){ set.seed(seed) }
  family = match.arg(family)
  co = corr^(abs(outer(c(1:p),c(1:p),"-"))); chco = chol(co)
  x.mat = matrix(rnorm(n*p),n,p)%*%t(chco)
  b.vec = seq(cf.max,cf.min,length.out=q)*(-1)^(1:q); b.vec = c(b.vec,rep(0,p-q))
  xb.vec = pmin(as.vector(x.mat%*%b.vec),700); exb.vec = exp(xb.vec)
  if(family=="gaussian"){ y.vec = xb.vec + rnorm(n) }
  if(family=="binomial"){ p.vec = exb.vec/(1+exb.vec); y.vec = rbinom(n,1,p.vec) }
  if(family=="poisson"){ m.vec = exb.vec; y.vec = rpois(n,m.vec) }
  if(family=="multinomial"){
    b.mat = matrix(0,p,k-1)
    for(i in 1:(k-1)) b.mat[1:q,i] = sample(seq(cf.max,cf.min,length.out=q)*(1)^(1:q))
    xb.mat = pmin(x.mat%*%b.mat,700); exb.mat = exp(xb.mat)
    p.mat = exb.mat / (1+rowSums(exb.mat))
    p.mat0 = cbind(p.mat,1-rowSums(p.mat))
    y.vec = rep(0,n)
    for(i in 1:n){ y.vec[i] = sample(1:k,1,prob=p.mat0[i,]) }
    b.vec = b.mat
  }
  if(family=="cox"){
    u = runif(n,0,1); y.vec = -log(u)/exb.vec
    cen = rbinom(n,prob=r,size=1); x.mat = cbind(x.mat,1-cen)
  }
  return(list(x.mat=x.mat,y.vec=y.vec,b.vec=b.vec,family=family))
}




