#' @exportPattern "^[[:alpha:]]+"
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @useDynLib ncpen
#'
#'
#' @title
#' non-convex penalized estimation for generalized linear models
#'
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
#' \dontrun{}
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
#' non-convex penalized estimation for generalized linear models
#'
#' @description
#' \code{ncpen} fits generalized linear models by using penalized maximum likelihood estimation.
#' It covers Various non-convex penalties such as SCAD and MCP for linear, logistic and poisson regression models.
#'
#' @param y.vec numeric vector; samples of dependent variable
#' @param x.mat numeric matrix; samples of independent variables
#' @param family character; model type; \code{gaussian} for linear,
#' \code{binomial} for logistic, and \code{poisson} for includes  ; defalut is "gaussian"
#' @param penalty character;

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
#' \dontrun{}
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
          cat("cv fold number:",fold,"\n")
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

coef.ncpen = function(ncpen.fit){
     num = dim(ncpen.fit$coefficients)[1]
     if(ncpen.fit$intercept==TRUE) rownames(ncpen.fit$coefficients) = c("intercept",rep("",num-1))
     return(ncpen.fit$coefficients)
}

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
#' plot.ncpen
#'
#' @description
#' \code{ncpen} fits generalized linear models by using penalized maximum likelihood estimation.
#' It covers Various non-convex penalties such as SCAD and MCP for linear, logistic and poisson regression models.
#'
#' @param y.vec numeric vector; samples of dependent variable
#' @param x.mat numeric matrix; samples of independent variables
#' @param family character; model type; \code{gaussian} for linear,
#' \code{binomial} for logistic, and \code{poisson} for includes  ; defalut is "gaussian"
#' @param penalty character;

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
#' \dontrun{}
#' @export
plot.ncpen = function(ncpen.fit,log.scale=FALSE,...){
     b.mat = ncpen.fit$coefficients[-1,]; lambda = ncpen.fit$lambda
     if(log.scale==TRUE) lambda = log(lambda)
     plot(lambda,b.mat[1,],type="n",ylim=c(min(b.mat),max(b.mat)),main="trace of coefficients",ylab="",...)
     for(pos in 1:dim(b.mat)[1]){ lines(lambda,b.mat[pos,],col=pos,lty=pos) }
}

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


coef.cv.ncpen = function(cv.ncpen.fit,type=c("error","deviance")){
     type = match.arg(type)
     if(type=="error"){
          return(cv.ncpen.fit$opt.ebeta)
     }
     if(type=="deviance"){
          return(cv.ncpen.fit$opt.dbeta)
     }
}

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


#-------------------------------------------------------------------------
