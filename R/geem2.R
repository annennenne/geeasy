#' Fit Generalized Estimating Equations
#' 
#' Calculate coefficients and nuisance parameters using generalized 
#' estimating equations.  Link and Variance functions can be 
#' specified by the user.  Similar to \code{\link{glm}}.
#' 
#' @param formula A formula expression similar to that for \code{\link{glm}}, 
#' of the form \code{response~predictors}.  An offset is allowed, as in \code{glm}.
#' 
#' @param id A vector identifying the clusters. By default, data are assumed 
#' to be sorted such that observations in a cluster are in consecutive rows 
#' and higher numbered rows in a cluster are assumed to be later.  
#' If NULL, then each observation is assigned its own cluster.
#' 
#' @param waves An integer vector identifying components of a cluster. 
#' For example, this could be a time ordering. If integers are skipped within 
#' a cluster, then dummy rows with weight 0 are added in an attempt to preserve
#' the correlation structure (except if \code{corstr = "exchangeable"} or 
#' \code{"independent"}). This can be skipped by setting \code{nodummy=TRUE}.
#' 
#' @param data An optional data frame containing the variables in the model.
#' 
#' @param family Will determine the link and variance functions.  The argument 
#' can be one of three options: a \code{family} object, a character string,
#'  or a list of functions. For more information on how to use \code{family} 
#'  objects, see details below. 
#'  
#' @param corstr A character string specifying the correlation structure.  
#'  Allowed structures are: \code{"independence"}, \code{"exchangeable"}, 
#'  \code{"ar1"}, \code{"m-dependent"}, \code{"unstructured"}, \code{"fixed"}, 
#'  and \code{"userdefined"}.  Any unique substring may be supplied.  
#'  If \code{"fixed"} or \code{"userdefined"}, then \code{corr.mat} must be 
#'  specified.  If \code{"m-dependent"}, then \code{Mv} is relevant.
#'  
#' @param Mv For \code{"m-dependent"}, the value for \code{m}. 
#'  
#' @param weights A vector of weights for each observation.  If an observation
#'   has weight 0, it is excluded from the calculations of any parameters.  
#'   Observations with a \code{NA} anywhere (even in variables not included
#'    in the model) will be assigned a weight of 0.  Note that these weights 
#'    are now the same as PROC GEE weights and not PROC GENMOD.
#'    
#' @param corr.mat The correlation matrix for \code{"fixed"}.  Matrix should
#'   be symmetric with dimensions >= the maximum cluster size.  If the correlation 
#'   structure is \code{"userdefined"}, then this is a matrix describing which 
#'   correlations are the same.
#'   
#' @param init.beta An optional vector with the initial values of beta.  If not 
#'  specified, then the intercept will be set to \code{InvLink(mean(response))}. 
#'   \code{init.beta} must be specified if not using an intercept.
#'   
#' @param init.alpha An optional scalar or vector giving the initial values for 
#'   the correlation.  If provided along with \code{Mv>1} or \code{unstructured} 
#'   correlation, then the user must ensure that the vector is of the appropriate 
#'   length.
#'   
#' @param init.phi An optional initial overdispersion parameter.  If not supplied, 
#'   initialized to 1.
#'   
#' @param scale.fix If set to \code{TRUE}, then the scale parameter is fixed at
#'    the value of \code{init.phi}.
#'    
#' @param nodummy If set to \code{TRUE}, then dummy rows will not be added
#'    based on the values in \code{waves}.
#'    
#' @param sandwich If \code{TRUE}, calculate robust variance.
#'    
#' @param useP If set to \code{FALSE}, do not use the n-p correction for 
#'    dispersion and correlation estimates, as in Liang and Zeger. This can be 
#'    useful when the number of observations is small, as subtracting p may yield 
#'    correlations greater than 1.
#'    
#' @param maxit Maximum number of iterations.
#'    
#' @param tol Tolerance in calculation of coefficients.
#' 
#' @details Users may specify functions for link and variance functions, but the
#'  functions must be vectorized functions.  See \code{\link{Vectorize}} for an easy
#'  way to vectorize functions.  \code{Vectorize} should be used sparingly, however, 
#'  as it can lead to fairly slow function calls.  Care must be taken to ensure
#'  that convergence is possible with non-standard functions.
#'  
#'  Offsets must be specified in the model formula, as in glm.
#'  
#'  For the \code{"userdefined"} correlation option, the function accepts a 
#'  matrix with consecutive integers.  \code{geem} only looks at the upper 
#'  triangle of the matrix.  Any entry given as 0 will be fixed at 0.  All
#'   entries given as 1 will be assumed to be the same as each other and will 
#'   be assumed to be possibly different from entries with a 2, and so on.
#'   
#'  If observations are dropped because they have a weight of 0, then the 
#'  denominator for the moment estimates of the correlation matrices are 
#'  calculated using the number of non-zero Pearson residuals for the 
#'  correlation structures \code{unstructured}, \code{userdefined} and 
#'  \code{m-dependent} with \code{Mv>1}.  Therefore residuals numerically 
#'  equal to 0 may cause problems in the calculation of correlation parameters.
#' 
#'  Concerning the \code{family} argument: If the supplied argument is a character 
#'  string, then the string should correspond to one of the family objects.
#'  In order to define a link function, a list must be created with the 
#'  components \code{(LinkFun, VarFun, InvLink, InvLinkDeriv)}, all of which are 
#'  vectorized functions.  If the components in the list are not named
#'   as \code{(LinkFun, VarFun, InvLink, InvLinkDeriv)}, then \code{geem} 
#'   assumes that the functions are given in that order.  LinkFun and VarFun 
#'   are the link and variance functions. InvLink and InvLinkDeriv are the inverse 
#'   of the link function and the derivative of the inverse of the link function 
#'   and so are decided by the choice of the link function.
#' 
#' @return An object of class \code{geem} representing the fit.
#' 
#' @author Lee McDaniel and Nick Henderson
#' 
#' @seealso \code{\link{glm}}, \code{\link{formula}}, \code{\link{family}}
#' 
#' @keywords models robust
#' 
#' @examples
#' 
#' ### Generated Negative Binomial Data
#' generatedata <- function(beta,alpha,gamma,X,T,n)  {
#'   mean.vec <- exp(crossprod(t(X),beta))
#'   y <- matrix(0,nrow=n,ncol=T)
#'   y[,1] <- rnbinom(n,mu = mean.vec[1],size=mean.vec[1]/gamma)
#'   for (i in 1:n)  {
#'       for (t in 2:T)  {
#'           innovation.mean <- mean.vec[t] - alpha*(sqrt(mean.vec[t]*mean.vec[t-1]))
#'           I <- rnbinom(1,mu= innovation.mean,size= innovation.mean/gamma)                              
#'           first.shape <- alpha*sqrt(mean.vec[t]*mean.vec[t-1])/gamma
#'           second.shape <- mean.vec[t-1]/gamma - first.shape
#'           u <- rbeta(1,shape1 = first.shape,shape2=second.shape)
#'           a <- rbinom(1,size=y[i,t-1],prob=u)
#'           y[i,t] = a + I
#'       }
#'   }
#'   longform <- c(t(y))
#'   print(apply(y,2,mean))
#'   simdata <- data.frame(count = longform, time = rep(X[,2],times=n),
#'                         subject=rep(c(1:n),each=T))
#'   return(simdata)
#'   }
#'   
#' X <- cbind(rep(1,5),c(-.5,-.25,0,.25,.5))
#' testdat <- generatedata(beta=c(1,.5),alpha=.2,gamma=.5,X=X,T=5,n=3000)
#' far1 <- geem(count~ time, id=subject ,data = testdat, family=poisson, 
#'              corstr="ar1")
#'              
#' ### Ohio respiratory data from geepack
#'  if(require(geepack)){
#'      data("ohio", package="geepack")
#'      resplogit <- geem(resp ~ age + smoke + age:smoke, id=id, data = ohio, 
#'                        family = binomial, corstr = "m-dep" , Mv = 1)
#'      LinkFun <- function(arg){qcauchy(arg)}
#'      InvLink <- function(arg){pcauchy(arg)}
#'      InvLinkDeriv <- function(arg){dcauchy(arg)}
#'      VarFun <- function(arg){arg*(1-arg)}
#'      FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)
#'      
#'      respcauchit <- geem(resp ~ age + smoke + age:smoke, id=id, data = ohio, 
#'                          family = FunList, corstr = "m-dep" , Mv=1)
#' }
#' 
#' ### Seizure data from geepack
#' if(require(geepack)){
#'     data("seizure", package="geepack")
#'     seiz.l <- reshape(seizure,
#'                       varying=list(c("base","y1", "y2", "y3", "y4")),
#'                       v.names="y", times=0:4, direction="long")
#'     seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
#'     seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
#'     seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)
#'     
#'     seiz <- geem(y~ x + trt + x:trt+ offset(log(t)), id=id,data = seiz.l, 
#'                  family = poisson, corstr = "exchangeable")
#' }
#' 
#' @export
geem2 <- function(formula, id, waves=NULL, data = parent.frame(),
                 family = gaussian, corstr = "independence", Mv = 1,
                 weights = NULL, corr.mat = NULL, #init.beta=NULL,
                # init.alpha=NULL, init.phi = 1, 
                 scale.fix = FALSE, nodummy = FALSE,
                 sandwich = TRUE, #useP = TRUE, #maxit = 20, #tol = 0.00001,
                 return = "geem",
                 control = geem.control()){
  
  ########################################################################
  #Check and prep input arguments ########################################
  ########################################################################
  
  call <- match.call()
  
  ### First, get all the relevant elements from the arguments
  dat <- model.frame(formula, data, na.action=na.pass)
  nn <- dim(dat)[1]
  
  #Make id / weight / waves argument so that they can match character
  #strings OR names provided in enclosing environment/data
  if(typeof(data) != "environment") {
    if(length(call$id) == 1){
      subj.col <- which(colnames(data) == call$id)
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      } else {
        id <- eval(call$id, envir=parent.frame())
      }
    } else if(is.null(call$id)) {
      id <- 1:nn
    }
    
    if(length(call$weights) == 1) {
      weights.col <- which(colnames(data) == call$weights)
      if(length(weights.col) > 0) {
        weights <- data[,weights.col]
      } else {
        weights <- eval(call$weights, envir=parent.frame())
      }
    } 
    
    if(length(call$waves) == 1) {
      waves.col <- which(colnames(data) == call$waves)
      if(length(waves.col) > 0) {
        waves <- data[,waves.col]
      } else {
        waves <- eval(call$waves, envir=parent.frame())
      }
    } else if(is.null(call$waves)) {
      waves <- NULL
    }
  }
  
  if(is.null(weights)) weights <- rep(1, nn)
  
  if(!is.integer(waves) & !is.null(waves)) stop("waves must be either an integer vector or NULL") 
  
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves

  
  #SORT THE DATA ACCORDING TO WAVES
  if(!is.null(waves)){
    dat <- dat[order(id, waves),]
  }else{
    dat <- dat[order(id),]
  }
  
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  # Note that we need to assign weight 0 to rows with NAs
  # in order to preserve the correlation structure
  na.inds <- NULL
  if(any(is.na(dat))){ #??? dat is constructed from formula + data => we only drop in NAs that we use. Fix documentation, it says otherwise!
    na.inds <- which(is.na(dat), arr.ind=T)
    #note: na.inds contains information on both row and columns of NAs - both are needed below. 
  }
  
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", 
               "fixed", "userdefined")
  #*# cor.match <- charmatch(corstr, cor.vec)
  corstr <- cor.vec[charmatch(corstr, cor.vec)]
  
  #!!!to do:
  # - add step to change corstr if settings simplify (for example "m-dependent" and Mv == 1 => "ar1")
  
  
  
  #!!! check what happens if length(cor.str) > 1 and check what happens if there is no match (old code below
  # NOT replaced yet):
  #*#  if(is.na(cor.match)){stop("Unsupported correlation structure")} #!!check if this still works
  #*#
  #*#   else if(cor.match == 0){
  #*#         stop("Ambiguous Correlation Structure Specification")
  #*#   }else{ 
  #*#   stop("Unsupported Correlation Structure")
  #*#   }
  
  if(!is.null(dat$waves)){
    wavespl <- split(dat$waves, dat$id)
    idspl <- split(dat$id, dat$id)
    
    maxwave <- rep(0, length(wavespl))
    incomp <- rep(0, length(wavespl))
    
    for(i in 1:length(wavespl)){
      maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 1
      if(maxwave[i] != length(wavespl[[i]])){
        incomp[i] <- 1
      }
    }
    
    #If there are gaps and correlation isn't independent or exchangeable
    #then we'll add some dummy rows
    if(!(corstr %in% c("independence", "exchangeable")) & 
       (sum(incomp) > 0) & !nodummy){
      #*# if( !is.element(cor.match,c(1,3)) & (sum(incomp) > 0) & !nodummy){
      dat <- dummyrows(formula, dat, incomp, maxwave, wavespl, idspl)
      id <- dat$id
      waves <- dat$waves
      weights <- dat$weights
    }
  }
  
  
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0 #!!! consider: should this be done in "weights" within dat instead?
    
    
    ##???? It looks like single imputation is swept in below - but maybe it's just filling in a placeholder
    ## value that won't be used? But why? those obs should be dropped anyway below?
    
    #*# for(i in unique(na.inds)[,2]){ #parenthesis mistake? gets same var multiple times
    for(i in unique(na.inds[,2])){
      if(is.factor(dat[,i])){
        dat[na.inds[,1], i] <- levels(dat[,i])[1]
      }else{
        dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
      }
    }
  }
  
  
  includedvec <- weights > 0
  
  
  inclsplit <- split(includedvec, id)
  

  
  #Find indexes to be dropped
  dropind <- NULL 
  
  if(corstr %in% c("independence", "exchangeable")) {
    dropind <- which(weights == 0)
  } else {
    #*# dropid <- NULL
    allobs <- TRUE
    if(any(!includedvec)){
      allobs <- FALSE
      
      uniqueid <- unique(id)
      n_uniqueid <- length(uniqueid)
      dodropid <- rep(FALSE, n_uniqueid)
      
      for(i in 1:n_uniqueid) {
        
        #Drop observation only if all obs from that individual (ID) are missing
        if(all(!inclsplit[[i]])){ 
          dodropid[i] <- TRUE
          #*# dropid <- c(dropid, unique(id)[i]) #!!!! DOESN'T WORK WHEN ID IS FACTOR - THEN IT IS COERCED TO NUM AND NOT MATCHED BELOW.
        }
      }
      #*#  if(length(dropid) > 0) {
      #*#    dropind <- which(is.element(id, dropid))
      #*#  }
      if (sum(dodropid) > 0) {
          dropind <- which(id %in% uniqueid[dodropid])  
      }
    }
  }
  if(length(dropind) > 0) { #drop indexes just found
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- id[-dropind]
  }

  #*#  nn <- dim(dat)[1] #note: now computed in geem.fit
 
  #find ud af hvor k er?
  
 #*# modterms <- terms(formula)
  
  X <- model.matrix(formula,dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)
  
  #add extra info to corstr if necessary
  corstr <- list(name = corstr, extra = NULL)
  if (corstr$name == "m-dependent") corstr$extra <- Mv
  if (corstr$name %in% c("fixed", "userdefined")) corstr$extra <- corr.mat
  
  
  ######################
  
  #call geem.fit here!
  #Args used:
  # - x = X
  # offset = offset
  # weights = weights
  # control = control
  # id = id
  # family = family
  # corstr = corstr 
  # allobs = allobs #!!! consider: can this be determined with the above information?
  
  
  #####################
 
  results <- geem.fit(x = X, y = Y, offset = offset, weights = weights,
                  control = control, id = id, family = family,
                  corstr = corstr, allobs = allobs, sandwich = sandwich)
  
  
  if (!results$converged) {
    warning("Did not converge")
  }
  if (results$unstable) {
    warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")
  }
  
  # }
  
  ##########################################################################################
  # Pack and return output #################################################################
  ##########################################################################################
  
  if (return == "geem") {
    # Create object of class geem with information about the fit
    dat <- model.frame(formula, data, na.action = na.pass) #!!! check - is this different than dat defined above?
    X <- model.matrix(formula, dat) #!!! check - is this different than X defined above?
 
    results$coefnames <- colnames(X)
    results$call <- call
    results$X <- X
    
#!!!!!!!!!!!!!!!!!!replace or leave out?    results$dropped <- dropid
    results$dropped <- "Put something here?"
     #NOTE: original geem sometimes included this slot, sometimes not, as it was
     #set to NULL if not used and hence (unintentionally?) dropped. 
     #Backwards compatibility may be in conflict with keeping a strict output structure. 
    
    results$terms <- terms(formula)
    results$y <- Y
    results$formula <- formula
    
    #reorder list to make it identical to previous structure
    old_geem_out_order <- c("beta", "phi", "alpha", "coefnames", 
                            "niter", "converged", "naiv.var", "var",
                            "call", "corr", "clusz", "FunList", 
                            "X", "offset", "eta", "weights", "terms",
                            "y", "biggest.R.alpha", "formula")
    results <-results[old_geem_out_order]
    
    class(results) <- "geem"
    
    #browser()
    
    return(results)
  } 
}


#####################################################################################
## Not exported below
#####################################################################################


### Simple moment estimator of dispersion parameter
#' @import Matrix
updatePhi <- function(YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW, useP){
  nn <- sum(includedlen)
  
  #  browser()
  
  
  resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu))
  
  phi <- (1/(sum(included)- useP * p))*crossprod(resid, resid)
  
  return(as.numeric(phi))
}

### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
updateBeta <- function(YY, XX, beta, offset, InvLinkDeriv, InvLink,
                       VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, W, included){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(XX%*%beta.new) + offset
    
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    hess <- crossprod( StdErr %*% dInvLinkdEta %*%XX, R.alpha.inv %*% W %*% StdErr %*%dInvLinkdEta %*% XX)
    esteq <- crossprod( StdErr %*%dInvLinkdEta %*%XX , R.alpha.inv %*% W %*% StdErr %*% as.matrix(YY - mu))
    
    #hess <- crossprod( StdErr %*% dInvLinkdEta %*%XX, included %*% R.alpha.inv  %*% W %*% StdErr %*%dInvLinkdEta %*% XX)
    #esteq <- crossprod( StdErr %*%dInvLinkdEta %*%XX , included %*% R.alpha.inv %*% W %*% StdErr %*% as.matrix(YY - mu))
    
    
    update <- solve(hess, esteq)
    
    
    beta.new <- beta.new + as.vector(update)
    
  }
  return(list(beta = beta.new, hess = hess))
}

### Calculate the sandiwch estimator as usual.
getSandwich <- function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv,
                        InvLink, VarFun, hessMat, StdErr, dInvLinkdEta,
                        BlockDiag, W, included){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX,  R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*%  StdErr %*% dInvLinkdEta %*% XX))
  #numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% XX, included %*% R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*% included %*% StdErr %*% dInvLinkdEta %*% XX))
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}



