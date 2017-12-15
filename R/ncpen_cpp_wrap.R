# @exportPattern "^[[:alpha:]]+"
#' @import RcppArmadillo
# @importFrom Rcpp evalCpp
#' @useDynLib ncpen
#'
#'
#'
#' @title
#' Fit a GLM with various nonconvex penalties
#'
#' @description
#' Fit a generalized linear model by penalized maximum likelihood estimation. The coefficients path is computed for the penalized regression model over a grid of values for the regularization parameter \eqn{\lambda}. Fit linear, logistic and poisson regression models with various non-convex penalties such as SCAD, MCP and clipped Lasso.
#'
#'
#' @param y.vec numeric vector; response vector
#' @param x.mat numeric matrix; design matrix; each row is an observation vector
#' @param family character; model type depending on the response; default is "gaussian"
#' @param penalty character; penalty function; default is "scad"
#' @param lambda numeric vector; a user-specified sequence of \code{lambda} values
#' @param n.lambda numeric; the number of lambda values; default is 100
#' @param r.lambda numeric; the smallest value for \code{lambda}, as a fraction of lambda.max (derived from data) for which all coefficients are zero
#' @param pen.weight numeric vector; penalty weights for each coefficient; a penalty weight can be zero in which the corresponding coefficient is always non-zero without shrinkage; note: the penalty weights are internally rescaled to sum to the number of variables, and the \code{lambda} sequence reflect this change
#' @param tau numeric; concavity parameter of the concave penalties (see reference); default is 3.7 for \code{scad}, 3 for \code{mcp}, 2 for \code{classo} and \code{sbridge}
#' @param gamma numeric; addtional tunning parameter for the \code{classo} and \code{sbridge}; default is 1e-6
#' @param ridge numeric; ridge effect (amount of ridge penalty); default is 1e-6
#' @param df.max numeric; upper bound for the number of nonzero coefficients
#' @param proj.min numeric; the minimum number of iterations which the projection is applied (see details); default is 100
#' @param iter.max numeric; maximum number of iterations; default is 1e+3
#' @param b.eps numeric; convergence threshold for L2 norms of coefficnets vector; default is 1e-6
#' @param k.eps numeric; convergence threshold for KKT conditions; default is 1e-4
#' @param x.standardize logical; standardization of \code{x.mat}, prior to fitting the model sequence. The coefficients are always returned on the original scale; default is \code{TRUE}
#' @param intercept logical; intercept in the model; default is \code{TRUE}
#'
#'
#' @details
#' The sequence of models indexed by the regularization parameter \code{lambda} is fit by the unified algorithm using concave convex procedure and coordiante descent algorithm. Note that the objective function is \deqn{ RSS / 2n +  penalty } for \code{family="gaussian"}, and \deqn{ negative log-likelihood / n + penalty } for \code{family="binomial"} or \code{family="poisson"}, where log-likelihood is computed with assuming the canonical link (logit for \code{binomial}; log for \code{poisson})
#'
#' The algorithm fits the coefficients in the active set using the projection method after \code{proj.min} iteration instead of cycling coordinates, which makes the algorithm fast and stable. (xxx to be modified)
#'
#'
#' @return An object with S3 class "ncpen" containing
#'   \item{warningings}{warnings from \code{native_cpp_ncpen_fun_} (xxx should be modified) }
#'   \item{family}{model type}
#'   \item{x.standardize}{standardization of \code{x.mat}}
#'   \item{intercept}{intercept in the model}
#'   \item{coefficients}{matrix of fitted coefficients for a \code{lambda} sequence; the number of rows is the number of coefficients (\code{ncol(x.mat)+1} if \code{intercept=T} and \code{ncol(x.mat)} if \code{intercept=F}) and the number of colums is equal to \code{nlambda}}
#'   \item{pen.weight}{penalty weights for each coefficient}
#'   \item{lambda}{actual sequence of \code{lambda} values used}
#'   \item{df}{the number of non-zero coefficients for each \code{lambda} value}
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
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' # 1. SCAD
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' coef(ncpen)
#' plot(ncpen)
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="regression")
#' gic.ncpen(ncpen,y.vec,x.mat)
#'
#' # 2. CLASSO
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian", penalty="classo")
#' plot(ncpen)
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="regression")
#'
#' # 3. TLP
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian", penalty="tlp")
#' plot(ncpen)
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="regression")
#'
#' ### Logistic regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="probability")
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="response")
#'
#' ### Poison regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="poisson")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="poisson")
#' predict(ncpen, new.x.mat=x.mat[1:20,],type="response")
#' gic.ncpen(ncpen,y.vec,x.mat)
#' plot(ncpen)
#' @export
ncpen = function(y.vec,x.mat,
                 family=c("gaussian","binomial","poisson"),
                 penalty=c("scad","mcp","lasso","classo","sridge","mbridge","mlog"),
                 lambda=NULL,n.lambda=1e+2,r.lambda=1e-3,w.lambda=NULL,
                 tau=5,gam=1e-6,ridge=1e-6,
                 df.max=50, proj.min=1e+2,
                 iter.max=1e+3,b.eps=1e-6,k.eps=1e-4,
                 x.standardize=TRUE,intercept=TRUE){
     family = match.arg(family)
     penalty = match.arg(penalty)
     n = dim(x.mat)[1]
     p = dim(x.mat)[2]
     if(is.null(lambda)){
          lambda = rep(-1,n.lambda)
     }
     if(is.null(w.lambda)){
          w.lambda = rep(1,p)
     } else {
          w.lambda = (p-sum(w.lambda==0))*w.lambda/sum(w.lambda)
     }
     if(sum(w.lambda==0)>n){
          stop("the number of zero (unpenalized) elements in w.lambda should be samller than the sample size")
     }
     if(sum(w.lambda==0)==p){
          stop("the number of zero (unpenalized) elements in w.lambda should be samller than the number of input variables")
     }
     ncpen.fit = native_cpp_ncpen_fun_(y.vec,
                                       x.mat,x.standardize,intercept,
                                       w.lambda,lambda, r.lambda,
                                       gam,tau,
                                       df.max,iter.max,b.eps,k.eps,proj.min,ridge,
                                       family, penalty);
     ret = list(warnings=as.vector(ncpen.fit$warnings),
                 family=family,x.standardize=x.standardize,intercept=intercept,
                 coefficients=ncpen.fit$b.mat,
                 w.lambda=w.lambda,                   ### w.lambda does not include the intercept
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
#'
#' @param y.vec numeric vector; the response vector as in \code{ncpen}
#' @param x.mat numeric matrix; the design matrix as in \code{ncpen}
#' @param n.fold numeric; the number of folds; default is 10; the smallest value allowable is 3
#' @param ... other arguments to \code{ncpen}
#'
#'
#' @details
#' The function runs \code{ncpen} \code{n.fold+1} times; the first to get the sequence of \code{lambda}, and then the remainder to compute the fit with each of the folds omitted. It provides the cross validated-error based on the squared-error loss and the deviance loss. (xxx the code may be modified)
#'
#'
#' @return An object with S3 class "cv.ncpen" containing
#'   \item{ncpen.fit}{the fitted \code{ncpen} object for the full data}
#'   \item{opt.ebeta}{the optimal coefficients vector selected by using the squared-error loss in the cross-valdation}
#'   \item{opt.dbeta}{the optimal coefficients vector selected by using the deviance loss in the cross-valdation}
#'   \item{cv.error}{the averaged cross-validated error for each value of \code{lambda}}
#'   \item{cv.deviance}{the averaged cross-validated deviance for each value of \code{lambda}}
#'   \item{elambda}{the actual \code{lambda} sequence used for computing cv error}
#'   \item{dlambda}{the actual \code{lambda} sequence used for computing cv deviance}
#'   \item{opt.elambda}{the optimal value of \code{lambda} based on cv error}
#'   \item{opt.dlambda}{the optimal value of \code{lambda} based on cv deviance}
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
#' set.seed(1234)
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian")
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
                    family="gaussian",penalty="scad",n.fold=10,
                    lambda=NULL,n.lambda=1e+2,r.lambda=1e-3,w.lambda=NULL,
                    tau=5,gam=1e-6,ridge=1e-6,
                    df.max=50, proj.min=1e+2,
                    iter.max=1e+3,b.eps=1e-6,k.eps=1e-4,
                    x.standardize=TRUE,intercept=TRUE){
     if(n.fold<2) stop("n.fold must be larger than 2")
     ncpen.fit = ncpen(y.vec,x.mat,family,penalty,lambda,n.lambda,r.lambda,w.lambda,
                       tau,gam,ridge,df.max,proj.min,iter.max,b.eps,k.eps,x.standardize,intercept)
     n = dim(x.mat)[1]
     p = dim(x.mat)[2]
     if(family=="gaussian"){
          cut = median(y.vec)
     } else {
          cut = 0
     }
     set1 = (1:n)[y.vec>median(y.vec)]
     set2 = (1:n)[y.vec<=median(y.vec)]
     f.list1 = suppressWarnings(split(sample(set1),1:n.fold))
     f.list2 = suppressWarnings(split(sample(set2),1:n.fold))
     err.list = list()
     dev.list = list()
     lambda = ncpen.fit$lambda
     w.lambda = ncpen.fit$w.lambda
     for(fold in 1:n.fold){
          # cat("cv fold number:",fold,"\n")
          tset = c(f.list1[[fold]],f.list2[[fold]])
          f.ncpen.fit = native_cpp_ncpen_fun_(y.vec[-tset],
                                              x.mat[-tset,],x.standardize,intercept,
                                              w.lambda,lambda, r.lambda,
                                              gam,tau,
                                              df.max,iter.max,b.eps,k.eps,proj.min,ridge,
                                              family, penalty);
          if(intercept==TRUE){ xb.mat = cbind(1,x.mat[tset,])%*%f.ncpen.fit$b.mat; xb.mat = pmin(xb.mat,700)
          } else { xb.mat = x.mat[tset,]%*%f.ncpen.fit$b.mat; xb.mat = pmin(xb.mat,700) }
          if(family=="gaussian"){ ny.mat = xb.mat }
          if(family=="binomial"){ ny.mat = ifelse(xb.mat>0,1,0); dev.list[[fold]] = 2*colSums(log(1+exp(xb.mat))-y.vec[tset]*xb.mat) }
          if(family=="poisson"){ ny.mat = exp(xb.mat); dev.list[[fold]] = 2*colSums(exp(xb.mat)-y.vec[tset]*xb.mat) }
          err.list[[fold]]= colSums((ny.mat-y.vec[tset])^2)
     }
     if(family=="gaussian") dev.list = err.list
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
                 opt.ebeta=ncpen.fit$coefficients[,err.opt],opt.dbeta=ncpen.fit$coefficients[,dev.opt],
                 cv.error=err.vec,cv.deviance=dev.vec,
                 elambda=lambda[1:eleng],dlambda=lambda[1:dleng],
                 opt.elambda=lambda[err.opt],opt.dlambda=lambda[dev.opt]
     ))
}

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
###### ncpen accessories  #####

#' @title
#' Extract the coefficients from a "ncpen" object
#'
#'
#' @description
#' This function returns the coefficients matrix for all lambda values
#'
#'
#' @param object Fitted "ncpen" model object
#'
#'
#' @return The coefficients matrix
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
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' coef(ncpen)
#'
#' @export
coef.ncpen = function(ncpen.fit){
     num = dim(ncpen.fit$coefficients)[1]
     if(ncpen.fit$intercept==TRUE) rownames(ncpen.fit$coefficients) = c("intercept",rep("",num-1))
     return(ncpen.fit$coefficients)
}

#' @title
#' Compute the GIC values for selection of the regularizatin parameter lambda
#'
#'
#' @description
#' This function provides the selection of the regularization parameter lambda based on the generlized information criterion(GIC) including AIC and BIC. It computes the GIC values at a grid of values for the regularization parametr lambda.
#'
#'
#' @param object Fitted "ncpen" model object
#' @param y.vec the response vector as in \code{ncpen}
#' @param x.mat the design matrix as in \code{ncpen}
#' @param df.weight the weight factor for various information critera; for examples, AIC if \code{df.weight=2}, BIC if \code{df.weight=log(n)}; default is BIC.
#' @param verbose logical; plot the GIC curve; default is \code{verbose=T}
#'
#'
#' @return The coefficients matrix
#'   \item{opt.beta }{the optimal coefficients vector selected by GIC }
#'   \item{lambda }{the sequence of lambda values in the "ncpen" object}
#'   \item{gic }{the GIC values at all lambda values}
#'   \item{opt.lambda }{the optimal lambda value}
#'   \item{plot }{the GIC curve}
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
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' gic.ncpen(ncpen,y.vec,x.mat,verbose=T)
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
#' Plot coefficients from a "ncpen" object
#'
#'
#' @description
#' Produces a plot of the coefficients paths for a fitted "ncpen" object
#'
#'
#' @param object Fitted "ncpen" model object
#' @param log.scale logical; log scale of horizontal axis(a sequence of lambda values); default is FALSE
#' @param ... other graphical parameters to \code{plot}
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
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5)
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' plot(ncpen,log.scale=F)
#'
#' @export
plot.ncpen = function(ncpen.fit,log.scale=FALSE,...){
     b.mat = ncpen.fit$coefficients[-1,]; lambda = ncpen.fit$lambda
     if(log.scale==TRUE) lambda = log(lambda)
     plot(lambda,b.mat[1,],type="n",ylim=c(min(b.mat),max(b.mat)),main="trace of coefficients",ylab="",...)
     for(pos in 1:dim(b.mat)[1]){ lines(lambda,b.mat[pos,],col=pos,lty=pos) }
}

#' @title
#' Make predictions from a "ncpen" object
#'
#'
#' @description
#' Similar to other predict methods, this function provides predictions from a fitted "ncpen" object.
#'
#'
#' @param object fitted "ncpen" model object
#' @param new.x.mat numeric matrix; matrix of new values at which predictions are to be made
#' @param type character; type of prediction required; "regression" returns the linear predictors; "probability" returns the fitted probabilities, which is only available for \code{family="binomial"}; "response" gives the fitted values for "gaussian", the fitted class thresholded by \code{cut} value for "binomial", and the fitted means for "poisson"
#' @param cut numeric; thresholding value of probability for logistic regression model; default is 0.5; this argument is not required for other models
#'
#'
#' @return the matrix of the fitted values depends on \code{type}, for all lambda values
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
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="gaussian")
#' predict(ncpen, new.x.mat=x.mat[1:20,], type="regression")
#'
#' ### Logistic regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' predict(ncpen, new.x.mat=x.mat[1:20,], type="probability")
#' predict(ncpen, new.x.mat=x.mat[1:20,], type="response")
#'
#' ### Poisson regression
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="poisson")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' ncpen = ncpen(y.vec=y.vec, x.mat=x.mat, family="poisson")
#' predict(ncpen, new.x.mat=x.mat[1:20,], type="regression")
#' predict(ncpen, new.x.mat=x.mat[1:20,], type="response")
#'
#' @export
predict.ncpen = function(ncpen.fit,new.x.mat=NULL,type=c("regression","probability","response"),cut=0.5){
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


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
###### ncpen.cv accessories  #####

#' @title
#' Make predictions from a "cv.ncpen" object
#'
#'
#' @description
#' This function makes predictions from a cross-validated ncpen model. It provides the optimal coefficients vector selected by the cross-validation method.
#'
#'
#' @param object fitted "cv.ncpen" model object
#' @param type character; type of cross-validated errors as in \code{cv.nepen}
#'
#'
#' @return the optimal coefficients vector selected by cross-validation method
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
#' \code{\link{cv.ncpen}}, \code{\link{cv.ncpen.plot}}
#'
#'
#' @examples
#'
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="binomial")
#' x.mat = s0$x.mat
#' y.vec = s0$y.vec
#'
#' cvfit = cv.ncpen(y.vec=y.vec, x.mat=x.mat, family="binomial")
#' coef.cv.ncpen(cvfit, type="deviance")
#'
#' @export
coef.cv.ncpen = function(cv.ncpen.fit,type=c("error","deviance")){
     type = match.arg(type)
     if(type=="error"){
          return(cv.ncpen.fit$opt.ebeta)
     }
     if(type=="deviance"){
          return(cv.ncpen.fit$opt.dbeta)
     }
}

#' @title
#' Plot cv curve from a "cv.ncpen" object
#'
#'
#' @description
#' Produces a plot of the cross-validated error curve from a fitted "cv.ncpen" object
#'
#'
#' @param object fitted "cv.ncpen" model object
#' @param type character; type of cross-validated errors as in \code{cv.nepen}
#' @param log.scale logical; log scale of horizontal axis(a sequence of lambda values); default is FALSE
#' @param ... other graphical parameters to \code{plot}
#'
#' @return the optimal coefficients vector selected by cross-validation method
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
#' @export
plot.cv.ncpen = function(cv.ncpen.fit,type=c("error","deviance"),log.scale=FALSE,...){
     type = match.arg(type)
     if(log.scale==TRUE){
          lambda = log(lambda)
     }
     if(type=="error"){
          plot(cv.ncpen.fit$elambda,cv.ncpen.fit$cv.error,main="trace of cv errors",ylab="",...)
     }
     if(type=="deviance"){
          plot(cv.ncpen.fit$dlambda,cv.ncpen.fit$cv.deviance,main="trace of cv deviances",ylab="",...)
     }

}

#' @title
#' Generate a dataset for simulations
#'
#' @description
#' Generate a synthetic dataset based on the correlation structure from generalized linear models.
#'
#'
#' @param n numeric; the number of samples
#' @param p numeric; the number of variables
#' @param q numeric; the number of nonzero coefficients
#' @param bmin numeric; value of the minimum coefficient
#' @param bmax numeric; value of the maximum coefficient
#' @param corr numeric; strength of correlations in the correlation structure
#' @param family character; model type; default is "gaussian"
#'
#'
#' @details
#' The design matrix for regression models is generated from the multivariate normal distribution with the correlation structure. Then the response variables are computed with a specific model based on the true coefficients. For details, see the reference.
#'
#'
#' @return An object with list class containing
#'   \item{x.mat }{\code{n} times \code{p} design matrix}
#'   \item{y.vec }{vector of responses}
#'   \item{b.vec }{vector of true coefficients}
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
#' set.seed(1234)
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian")
#' head(s0$x.mat)
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' s0 = sam.gen.fun(n=100,p=20,q=10,bmin=0.2,bmax=0.5,corr=0.7,family="binomial")
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' s0 = sam.gen.fun(n=100,p=20,q=5,bmin=0.5,bmax=1,corr=0.3,family="poisson")
#' head(s0$y.vec)
#' head(s0$b.vec)
#'
#' @export
sam.gen.fun = function(n=100,p=50,q=10,bmin=0.5,bmax=1,corr=0.5,family="gaussian"){
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
