#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline lines plot
#' @importFrom stats formula median model.matrix rbinom rnorm rpois
#' @useDynLib ncpen, .registration = TRUE
#'
#'
#'
#' @title
#' Fits a generalized linear model (GLM) with various nonconvex penalties
#'
#' @description
#' Fits a generalized linear model by penalized maximum likelihood estimation.
#' The coefficients path is computed for the penalized regression model over a grid of values for the regularization parameter \eqn{\lambda}.
#' Fits gaussian (linear), binomial (logistic) and poisson regression models with various non-convex penalties such as SCAD, MCP and clipped Lasso.
#'
#'
#' @param y.vec (numeric vector) response vector.
#' @param x.mat (numeric matrix) design matrix. Each row is an observation vector.
#' @param family (character) regression model. Default is "\code{gaussian}".
#' @param penalty (character) penalty function. Default is "\code{scad}".
#' @param lambda (numeric vector): user-specified sequence of \code{lambda} values.
#' @param n.lambda (numeric) the number of \code{lambda} values. Default is 100.
#' @param r.lambda (numeric) ratio of the smallest value for \code{lambda} to \code{lambda.max} (which derived from data) for which all coefficients are zero. Default is 1e-3.
#' @param pen.weight (numeric vector) penalty weights for each coefficient. If a penalty weight is set to zero,
#' the corresponding coefficient is always non-zero without shrinkage.
#' Note: the penalty weights are internally rescaled to sum to the number of variables, and the \code{lambda} sequence reflects this change.
#' @param tau (numeric) concavity parameter of the concave penalties (see reference). Default is 3.7 for \code{scad}, 3 for \code{mcp}, 2 for \code{classo} and \code{sridge}, 0.1 for \code{tlp}, \code{mbridge} and \code{mlog}.
#' @param gamma (numeric) additional tuning parameter for the \code{classo} and \code{sbridge}. Default value is 1e-6.
#' @param ridge (numeric) ridge effect (amount of ridge penalty). Default value is 1e-6.
#' @param df.max (numeric) the maximum number of nonzero coefficients. Default is 50.
#' @param proj.min (numeric) the minimum number of iterations which will be applied to projections (see details). Default value is 50.
#' @param iter.max (numeric) maximum number of iterations. Default valu eis 1e+3.
#' @param b.eps (numeric) convergence threshold for \eqn{L2} norms of coefficients vector. Default value is 1e-7.
#' @param k.eps (numeric) convergence threshold for KKT conditions. Default value is 1e-6.
#' @param x.standardize (logical) whether to standardize the \code{x.mat} prior to fitting the model.
#' The estimated coefficients are always restored to the original scale. Default value is \code{TRUE}.
#' @param intercept (logical) whether to include an intercept in the model. Default value is \code{TRUE}.
#'
#'
#' @details
#' The sequence of models indexed by the regularization parameter \code{lambda} is fit by the unified algorithm
#' using concave convex procedure and coordinate descent algorithm.
#' Note that the objective function is \deqn{ RSS / 2n +  penalty } for \code{family="gaussian"},
#' and \deqn{(negative  log-likelihood) / n + penalty }
#' for \code{family="binomial"} or \code{family="poisson"}, where log-likelihood is computed with assuming the canonical link
#' (logit for \code{binomial}; log for \code{poisson}).
#'
#' The algorithm fits the coefficients in the active set using the projection method
#' after \code{proj.min} iteration instead of cycling coordinates, which makes the algorithm fast and stable.
#'
#'
#' @return An object with S3 class \code{ncpen}.
#'   \item{family}{regression model.}
#'   \item{x.standardize}{flag for standardization of \code{x.mat}.}
#'   \item{intercept}{flag for an intercept in the model.}
#'   \item{coefficients}{a matrix of fitted coefficients for a \code{lambda} sequence.
#'   The number of rows is same as the number of coefficients (\code{ncol(x.mat)+1} if \code{intercept=TRUE} and
#'   \code{ncol(x.mat)} if \code{intercept=FALSE}). The number of columns is equal to \code{nlambda}.}
#'   \item{pen.weight}{penalty weights for each coefficient.}
#'   \item{lambda}{sequence of \code{lambda} values used.}
#'   \item{df}{the number of non-zero coefficients for each \code{lambda} value.}
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{plot.ncpen}}, \code{\link{coef.ncpen}}, \code{\link{cv.ncpen}}
#'
#'
#' @examples
#'
#' ### Linear regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' # 1. SCAD
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' coef(fit)
#' plot(fit)
#' predict(fit, new.x.mat=x.mat[1:20,],type="regression")
#' gic.ncpen(fit,y.vec,x.mat)
#'
#' # 2. CLASSO
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian", penalty="classo")
#' plot(fit)
#' predict(fit, new.x.mat=x.mat[1:20,],type="regression")
#'
#' # 3. TLP
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian", penalty="tlp")
#' plot(fit)
#' predict(fit, new.x.mat=x.mat[1:20,],type="regression")
#'
#' ### Logistic regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' predict(fit, new.x.mat=x.mat[1:20,],type="probability")
#' predict(fit, new.x.mat=x.mat[1:20,],type="response")
#'
#' ### Poison regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="poisson", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="poisson")
#' predict(fit, new.x.mat=x.mat[1:20,],type="response")
#' gic.ncpen(fit,y.vec,x.mat)
#' plot(fit)
#' @export
ncpen = function(y.vec,x.mat,
                 family=c("gaussian","binomial","poisson"),
                 penalty=c("scad","mcp","tlp","lasso","classo","sridge","mbridge","mlog"),
                 lambda=NULL,n.lambda=1e+2,r.lambda=1e-3,
                 pen.weight=NULL,
                 tau=switch(penalty,scad=3.7,mcp=3,tlp=0.1,lasso=1,classo=2,sridge=2,mbridge=0.1,mlog=0.1),
                 gamma=1e-6,ridge=1e-6,
                 df.max=50, proj.min=50,
                 iter.max=1e+3,b.eps=1e-7,k.eps=1e-6,
                 x.standardize=TRUE,intercept=TRUE){
     family = match.arg(family)
     penalty = match.arg(penalty)
     n = dim(x.mat)[1]
     p = dim(x.mat)[2]
     if(is.null(lambda)){
          lambda = rep(-1,n.lambda)
     }
     if(is.null(pen.weight)){
          pen.weight = rep(1,p)
     } else {
          pen.weight = (p-sum(pen.weight==0))*pen.weight/sum(pen.weight)
     }
     if(sum(pen.weight==0)>n){
          stop("the number of zero (unpenalized) elements in pen.weight should be samller than the sample size")
     }
     if(sum(pen.weight==0)==p){
          stop("the number of zero (unpenalized) elements in pen.weight should be samller than the number of input variables")
     }
     ncpen.fit = native_cpp_ncpen_fun_(y.vec,
                                       x.mat,x.standardize,intercept,
                                       pen.weight,lambda,r.lambda,
                                       gamma,tau,
                                       df.max,iter.max,b.eps,k.eps,proj.min,ridge,
                                       family, penalty);

     # Throw warning messages from cpp
     if(ncpen.fit$warnings[1] == 1) {
          warning("(ncpen: warning code 1) increase p.max.");
     }

     if(ncpen.fit$warnings[2] == 1) {
          warning("(ncpen: warning code 2) increase r.eff.");
     }


     ret = list(family=family,x.standardize=x.standardize,intercept=intercept,
                coefficients=ncpen.fit$b.mat,
                pen.weight=pen.weight,                   ### pen.weight does not include the intercept
                lambda=as.vector(ncpen.fit$lam.vec),
                df=as.vector(ncpen.fit$d.vec));

     class(ret) = "ncpen";
     return(ret);
}

#' @title
#' Cross-validation for ncpen
#'
#'
#' @description
#' Performs k-fold cross-validation for noncovex penalized regression models over a sequence of the regularization parameter lambda.
#'
#' @param y.vec (numeric vector) response vector.
#' @param x.mat (numeric matrix) design matrix. Each row is an observation vector.
#' @param family (character) regression model. Default is "\code{gaussian}".
#' @param penalty (character) penalty function. Default is "\code{scad}".
#' @param lambda (numeric vector): user-specified sequence of \code{lambda} values.
#' @param n.lambda (numeric) the number of \code{lambda} values. Default is 100.
#' @param r.lambda (numeric) ratio of the smallest value for \code{lambda} to \code{lambda.max} (which derived from data) for which all coefficients are zero. Default is 1e-3.
#' @param pen.weight (numeric vector) penalty weights for each coefficient. If a penalty weigth is set to zero,
#' the corresponding coefficient is always non-zero without shrinkage.
#' Note: the penalty weights are internally rescaled to sum to the number of variables, and the \code{lambda} sequence reflects this change.
#' @param tau (numeric) concavity parameter of the concave penalties (see reference). Default is 3.7 for \code{scad}, 3 for \code{mcp}, 2 for \code{classo} and \code{sridge}, 0.1 for \code{tlp}, \code{mbridge} and \code{mlog}.
#' @param gamma (numeric) addtional tuning parameter for the \code{classo} and \code{sbridge}. Default value is 1e-6.
#' @param ridge (numeric) ridge effect (amount of ridge penalty). Default value is 1e-6.
#' @param df.max (numeric) the maximum number of nonzero coefficients. Default is 50.
#' @param proj.min (numeric) the minimum number of iterations which will be applied to projections (see details). Default value is 50.
#' @param iter.max (numeric) maximum number of iterations. Default value is 1e+3.
#' @param b.eps (numeric) convergence threshold for \eqn{L2} norms of coefficients vector. Default value is 1e-7.
#' @param k.eps (numeric) convergence threshold for KKT conditions. Default value is 1e-6.
#' @param x.standardize (logical) whether to standardize the \code{x.mat} prior to fitting the model.
#' The estimated coefficients are always restored to the original scale. Default value is \code{TRUE}.
#' @param intercept (logical) whether to include an intercept in the model. Default value is \code{TRUE}.
#' @param n.fold (numeric) the number of folds. Default value is 10. It should be 3 or greater.
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
#'   \item{opt.dbeta}{the optimal coefficients vector selected by using the deviance loss in the cross-valdation.}
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
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{ncpen}}, \code{\link{plot.cv.ncpen}}, \code{\link{coef.cv.ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' cvfit = cv.ncpen(y.vec=y.vec,x.mat=x.mat,family="gaussian",n.fold=10)  # not run !!!
#' coef.cv.ncpen(cvfit)
#' plot.cv.ncpen(cvfit)
#' fit = cvfit$ncpen.fit
#' opt = which(cvfit$opt.elambda==fit$lambda)
#' coef(fit)[,opt]
#'
#' @export
cv.ncpen = function(y.vec,x.mat,
                    family=c("gaussian","binomial","poisson"),
                    penalty=c("scad","mcp","tlp","lasso","classo","sridge","mbridge","mlog"),
                    n.fold=10,
                    lambda=NULL,n.lambda=1e+2,r.lambda=1e-3,pen.weight=NULL,
                    tau=switch(penalty,scad=3.7,mcp=3,tlp=0.1,lasso=1,classo=2,sridge=2,mbridge=0.1,mlog=0.1),
                    gamma=1e-6,ridge=1e-6,
                    df.max=50, proj.min=50,
                    iter.max=1e+3,b.eps=1e-7,k.eps=1e-6,
                    x.standardize=TRUE,intercept=TRUE){
     family = match.arg(family)
     penalty = match.arg(penalty)
     if(n.fold<2) stop("n.fold must be larger than 2")
     ncpen.fit = ncpen(y.vec,x.mat,family,penalty,lambda,n.lambda,r.lambda,pen.weight,
                       tau,gamma,ridge,df.max,proj.min,iter.max,b.eps,k.eps,x.standardize,intercept)
     n = dim(x.mat)[1]
     p = dim(x.mat)[2]
     if(family=="gaussian"){ cut = median(y.vec) } else { cut = 0 }
     set1 = (1:n)[y.vec>median(y.vec)]
     set2 = (1:n)[y.vec<=median(y.vec)]
     f.list1 = suppressWarnings(split(sample(set1),1:n.fold))
     f.list2 = suppressWarnings(split(sample(set2),1:n.fold))
     err.list = list()
     dev.list = list()
     lambda = ncpen.fit$lambda
     pen.weight = ncpen.fit$pen.weight
     for(fold in 1:n.fold){
          cat("cv fold number:",fold,"\n")
          tset = c(f.list1[[fold]],f.list2[[fold]])
          f.ncpen.fit = native_cpp_ncpen_fun_(y.vec[-tset],
                                              x.mat[-tset,],x.standardize,intercept,
                                              pen.weight,lambda,r.lambda,
                                              gamma,tau,
                                              df.max,iter.max,b.eps,k.eps,proj.min,ridge,
                                              family, penalty);
          if(intercept==TRUE){
               xb.mat = cbind(1,x.mat[tset,])%*%f.ncpen.fit$b.mat;
               xb.mat = pmin(xb.mat,700)
          } else {
               xb.mat = x.mat[tset,]%*%f.ncpen.fit$b.mat;
               xb.mat = pmin(xb.mat,700)
          }
          if(family=="gaussian"){
               ny.mat = xb.mat
          }
          if(family=="binomial"){
               ny.mat = ifelse(xb.mat>0,1,0);
               dev.list[[fold]] = 2*colSums(log(1+exp(xb.mat))-y.vec[tset]*xb.mat)
          }
          if(family=="poisson"){
               ny.mat = exp(xb.mat);
               dev.list[[fold]] = 2*colSums(exp(xb.mat)-y.vec[tset]*xb.mat)
          }
          err.list[[fold]]= colSums((ny.mat-y.vec[tset])^2)
     }
     if(family=="gaussian"){
          dev.list = err.list
     }
     ### response prediction
     eleng = min(sapply(err.list,length))
     err.mat = matrix(0,n.fold,eleng)
     for(fold in 1:n.fold){ err.mat[fold,] = err.list[[fold]][1:eleng] }
     err.vec = colMeans(err.mat); err.opt = which.min(err.vec)
     ### deviance prediction
     dleng = min(sapply(dev.list,length))
     dev.mat = matrix(0,n.fold,dleng)
     for(fold in 1:n.fold){ dev.mat[fold,] = dev.list[[fold]][1:dleng] }
     dev.vec = colMeans(dev.mat); dev.opt = which.min(dev.vec)

     return(list(ncpen.fit=ncpen.fit,
                 opt.ebeta=ncpen.fit$coefficients[,err.opt],
                 opt.dbeta=ncpen.fit$coefficients[,dev.opt],
                 cv.error=err.vec,
                 cv.deviance=dev.vec,
                 elambda=lambda[1:eleng],
                 dlambda=lambda[1:dleng],
                 opt.elambda=lambda[err.opt],
                 opt.dlambda=lambda[dev.opt]
     ))
}


################################################################################################################
################################################################################################################
#######   ncpen accessories  #######

#' @title
#' Extract the coefficients from an \code{ncpen} object
#'
#'
#' @description
#' This function returns the coefficients matrix for all lambda values.
#'
#'
#' @param object Fitted \code{ncpen} object.
#' @param ... Other parameters to coef. Not supported.
#'
#'
#' @return The coefficients \code{\link{matrix}}.
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{ncpen}}, \code{\link{plot.ncpen}}, \code{\link{predict.ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5, seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' coef(fit)
#'
#' @export coef.ncpen
#' @export
coef.ncpen = function(object, ...){
     # S3 naming standard
     ncpen.fit = object;

     num = dim(ncpen.fit$coefficients)[1]
     if(ncpen.fit$intercept==TRUE) rownames(ncpen.fit$coefficients) = c("intercept",rep("",num-1))
     return(ncpen.fit$coefficients)
}


#' @title
#' Compute the GIC values for the selection of the regularization parameter lambda.
#'
#'
#' @description
#' This function provides the selection of the regularization parameter lambda based
#' on the generalized information criterion (GIC) including AIC and BIC.
#' It computes the GIC values at a grid of values for the regularization parameter lambda.
#'
#' @param ncpen.fit Fitted \code{ncpen} model object.
#' @param y.vec the response vector.
#' @param x.mat the design matrix.
#' @param df.weight the weight factor for various information criteria. For example, AIC if \code{df.weight=2},
#' BIC if \code{df.weight=log(n)}. Default is BIC.
#' @param verbose (logical) whether to plot the GIC curve. Default is \code{verbose=TRUE}.
#'
#'
#' @return The coefficients \code{\link{matrix}}.
#'   \item{opt.beta}{the optimal coefficients \code{\link{vector}} selected by GIC.}
#'   \item{lambda}{the sequence of lambda values in the \code{ncpen} object.}
#'   \item{gic}{the GIC values for all lambda values.}
#'   \item{opt.lambda}{the optimal lambda value.}
#'   \item{plot}{the GIC curve.}
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#' Kim, Y., Kwon, S. and Choi, H. (2012). Consistent Model Selection Criteria on High Dimensions. \emph{Journal of Machine Learning Research}, \bold{13}, 1037-1057.
#'
#'
#' @seealso
#' \code{\link{ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5, seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' gic.ncpen(fit,y.vec,x.mat,verbose=TRUE)
#'
#' @export
gic.ncpen = function(ncpen.fit,y.vec,x.mat,df.weight=log(length(y.vec)),verbose=TRUE){
     n = length(y.vec);
     if(ncpen.fit$intercept==FALSE){
          xb.mat = x.mat%*%ncpen.fit$coefficients
     } else {
          xb.mat = cbind(1,x.mat)%*%ncpen.fit$coefficients
     }
     if(ncpen.fit$family=="gaussian"){ dev.vec = n*log(colSums((y.vec-xb.mat)^2)/n) }
     if(ncpen.fit$family=="binomial"){ dev.vec = 2*colSums(log(1+exp(xb.mat))-y.vec*xb.mat)  }
     if(ncpen.fit$family=="poisson") { dev.vec = 2*colSums(exp(xb.mat)-y.vec*xb.mat) }
     df.vec = colSums(ncpen.fit$coefficients!=0)
     gic.vec = dev.vec + df.weight*df.vec
     if(verbose==TRUE){
          plot(ncpen.fit$lambda,gic.vec,xlab="lambda",ylab="",main="information criterion curve",type="b")
     }
     opt = which.min(gic.vec)
     return(list(opt.beta=ncpen.fit$coefficients[,opt],lambda=ncpen.fit$lambda,gic=gic.vec,opt.lambda=ncpen.fit$lambda[opt]))
}

#' @title
#' Plots coefficients from an \code{ncpen} object.
#'
#'
#' @description
#' Produces a plot of the coefficients paths for a fitted \code{ncpen} object.
#'
#'
#' @param x Fitted \code{ncpen} model object.
#' @param log.scale (logical) log scale of horizontal axis (a sequence of lambda values). Default value is FALSE.
#' @param ... other graphical parameters to \code{\link{plot}}
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5, seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' plot(fit,log.scale=FALSE)
#'
#' @export plot.ncpen
#' @export
plot.ncpen = function(x,log.scale=FALSE,...){
     # S3 standard naming
     ncpen.fit = x;

     b.mat = ncpen.fit$coefficients[-1,]; lambda = ncpen.fit$lambda
     if(log.scale==TRUE) lambda = log(lambda)
     plot(lambda,b.mat[1,],type="n",ylim=c(min(b.mat),max(b.mat)),main="trace of coefficients",ylab="",...)
     for(pos in 1:dim(b.mat)[1]){ lines(lambda,b.mat[pos,],col=pos,lty=pos) }
}

#' @title
#' Make predictions from an \code{ncpen} object.
#'
#'
#' @description
#' This function provides predictions from a fitted \code{ncpen} object.
#'
#'
#' @param object fitted \code{ncpen} object.
#' @param new.x.mat (numeric \code{\link{matrix}}). A matrix of new observations at which predictions are to be made.
#' @param type (character) type of prediction.
#' \code{"regression"} returns the linear predictors;
#' \code{"probability"} returns the fitted probabilities which is only available for \code{family="binomial"};
#' \code{"response"} returns followings depending on the models: the fitted values for \code{"gaussian"},
#' fitted class using \code{cut} value for \code{"binomial"}, and fitted means for \code{"poisson"}.
#' @param cut (numeric) threshold value of probability for logistic regression model. Default value is 0.5.
#' This argument is only required for logistic regression (binomial).
#' @param ... Other parameters to prediction. Not supported.
#'
#'
#' @return the \code{\link{matrix}} of the fitted values depending on \code{type} for all lambda values.
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{ncpen}}
#'
#'
#' @examples
#'
#' ### Linear regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' predict(fit, new.x.mat=x.mat[1:20,], type="regression")
#'
#' ### Logistic regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' predict(fit, new.x.mat=x.mat[1:20,], type="probability")
#' predict(fit, new.x.mat=x.mat[1:20,], type="response")
#'
#' ### Poisson regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="poisson")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' fit = ncpen(y.vec=y.vec, x.mat=x.mat, family="poisson")
#' predict(fit, new.x.mat=x.mat[1:20,], type="regression")
#' predict(fit, new.x.mat=x.mat[1:20,], type="response")
#'
#' @export predict.ncpen
#' @export
predict.ncpen = function(object,new.x.mat=NULL,type=c("regression","probability","response"),cut=0.5, ...){
     # S3 standard naming
     ncpen.fit = object;

     if(is.null(new.x.mat)){
          stop("'new.x.mat' option should be supplied for prediction")
     }
     if(ncpen.fit$intercept==TRUE){
          xb.mat = cbind(1,new.x.mat)%*%ncpen.fit$coefficients;
          xb.mat = pmin(xb.mat,700)
     } else {
          xb.mat = new.x.mat%*%ncpen.fit$coefficients;
          xb.mat = pmin(xb.mat,700)
     }
     type = match.arg(type)
     if(type=="regression"){
          return(pred.mat=xb.mat)
     }
     if(type=="probability"){
          if(ncpen.fit$family!="binomial"){
               stop("'family' option in 'ncpen' should be 'binomial'")
          } else {
               exb.mat = exp(xb.mat)
               p.mat = exb.mat/(1+exb.mat)
          }
          return(pred.mat=p.mat)
     }
     if(type=="response"){
          if(ncpen.fit$family=="gaussian"){
               ny.mat = xb.mat
          }
          if(ncpen.fit$family=="binomial") {
               exb.mat = exp(xb.mat)
               p.mat = exb.mat/(1+exb.mat)
               ny.mat = 1*(p.mat>cut)
          }
          if(ncpen.fit$family=="poisson") {
               ny.mat = exp(xb.mat)
          }
          return(pred.mat=ny.mat)
     }
}


################################################################################################################
################################################################################################################
###### ncpen.cv accessories  #######

#' @title
#' Extracts the optimal vector of coefficients from a \code{cv.ncpen} object.
#'
#'
#' @description
#' This function returns the optimal vector of coefficients.
#'
#'
#' @param object fitted \code{cv.ncpen} object.
#' @param type (character) a cross-validated error type which is either "error" or "deviance". Each error type is defined in \code{\link{cv.ncpen}}.
#' @param ... Other arguments to coef. Not supported.
#'
#'
#' @return the optimal coefficients vector selected by cross-validation method.
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{cv.ncpen}}, \code{\link{plot.cv.ncpen}}
#'
#'
#' @examples
#'
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial", seed = 1234)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' cvfit = cv.ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' coef.cv.ncpen(cvfit, type="deviance")
#'
#' @export coef.cv.ncpen
#' @export
coef.cv.ncpen = function(object,type=c("error","deviance"), ...){
     # S3 naming standard
     cvfit = object;

     type = match.arg(type)
     if(type=="error"){
          return(cvfit$opt.ebeta)
     }
     if(type=="deviance"){
          return(cvfit$opt.dbeta)
     }
}

#' @title
#' Plot cv curve from a \code{cv.ncpen} object
#'
#'
#' @description
#' Produces a plot of the cross-validated error curve from a fitted \code{cv.ncpen} object.
#'
#'
#' @param x fitted \code{cv.ncpen} object.
#' @param type (character) a cross-validated error type which is either "error" or "deviance". Each error type is defined in \code{\link{cv.ncpen}}.
#' @param log.scale (logical) log scale of horizontal axis (a sequence of lambda values). Default value is FALSE.
#' @param ... other graphical parameters to \code{plot}.
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{cv.ncpen}}, \code{\link{coef.cv.ncpen}}
#'
#'
#' @examples
#'
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' cvfit = cv.ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' plot.cv.ncpen(cvfit, type="deviance")
#'
#' @export plot.cv.ncpen
#' @export
plot.cv.ncpen = function(x,type=c("error","deviance"),log.scale=FALSE, ...){
     # For S3 generic/method consistency
     cvfit = x;

     type = match.arg(type)
     if(log.scale==TRUE){
          lambda = log(lambda)
     }
     if(type=="error"){
          plot(cvfit$elambda,cvfit$cv.error,main="trace of cv errors",ylab="",xlab="lambda",...)
          abline(v=cvfit$opt.elambda,col="gray50")
     }
     if(type=="deviance"){
          plot(cvfit$dlambda,cvfit$cv.deviance,main="trace of cv deviances",ylab="",xlab="lambda",...)
          abline(v=cvfit$opt.dlambda,col="gray50")
     }
}


#' @title
#' Generate a simulated dataset.
#'
#' @description
#' Generate a synthetic dataset based on the correlation structure from generalized linear models.
#'
#'
#' @param n (numeric) the number of samples.
#' @param p (numeric) the number of variables.
#' @param q (numeric) the number of nonzero coefficients.
#' @param bmin (numeric) value of the minimum coefficient.
#' @param bmax (numeric) value of the maximum coefficient.
#' @param corr (numeric) strength of correlations in the correlation structure.
#' @param family (character) model type. Default is \code{"gaussian"}.
#' @param seed (numeric) seed number for random generation. If set to NA, no seed will be applied. Default value is NA.
#'
#'
#' @details
#' A design matrix for regression models is generated from the multivariate normal distribution with the correlation structure.
#' Then the response variables are computed with a specific model based on the true coefficients. For details, see the reference.
#'
#'
#' @return An object with list class containing
#'   \item{x.mat}{  \code{n} times \code{p} design matrix.}
#'   \item{y.vec}{  vector of responses.}
#'   \item{b.vec}{  vector of true coefficients.}
#'
#'
#' @author Dongshin Kim, Sunghoon Kwon, Sangin Lee
#'
#' Maintainer: Dongshin Kim<dongshin.kim@outlook.com>, Sunghoon Kwon<shkwon0522@gmail.com>
#'
#'
#' @references
#' Kim, D., Kwon, S. and Lee, S. (2017). A unified algorithm for various penalized regression models: \bold{R} Package \bold{ncpen}.
#'
#'
#' @seealso
#' \code{\link{ncpen}}
#'
#'
#' @examples
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian", seed = 1234)
#' head(s0$x.mat)
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.2,bmax=0.5,corr=0.7,family="binomial", seed = 1234)
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' s0 = sam.gen.fun(n=100,p=20,q=5,bmin=0.5,bmax=1,corr=0.3,family="poisson", seed = 1234)
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' @export
sam.gen.fun = function(n=100,p=50,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian", seed = NA){
     if( !(is.na(seed) | is.null(seed)) ) {
          set.seed(seed);
     }

     co = corr^(abs(outer(c(1:p),c(1:p),"-"))); chco = chol(co)
     x.mat = matrix(rnorm(n*p),n,p)%*%t(chco)
     b.vec = seq(bmax,bmin,length.out=q)*(-1)^(1:q); b.vec = c(b.vec,rep(0,p-q))
     xb.vec = as.vector(x.mat%*%b.vec); exb.vec = pmin(exp(xb.vec),exp(700))
     if(family=="gaussian"){ y.vec = xb.vec + rnorm(n) }
     if(family=="binomial"){ p.vec = exb.vec/(1+exb.vec); y.vec = rbinom(n,1,p.vec) }
     if(family=="poisson"){ m.vec = exb.vec; y.vec = rpois(n,m.vec) }
     return(list(x.mat=x.mat,y.vec=y.vec,b.vec=b.vec))
}

#-------------------------------------------------------------------------
