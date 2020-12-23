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
geem <- function(formula, id, waves=NULL, data = parent.frame(),
                 family = gaussian, corstr = "independence", Mv = 1,
                 weights = NULL, corr.mat = NULL, init.beta=NULL,
                 init.alpha=NULL, init.phi = 1, scale.fix = FALSE, nodummy = FALSE,
                 sandwich = TRUE, useP = TRUE, maxit = 20, tol = 0.00001){
  call <- match.call()
  
  famret <- getfam(family)
  
  if(inherits(famret, "family")){
    LinkFun <- famret$linkfun
    InvLink <- famret$linkinv
    VarFun <- famret$variance
    InvLinkDeriv <- famret$mu.eta
  }else{
    LinkFun <- famret$LinkFun
    VarFun <- famret$VarFun
    InvLink <- famret$InvLink
    InvLinkDeriv <- famret$InvLinkDeriv
  }
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  useP <- as.numeric(useP)
  
  ### First, get all the relevant elements from the arguments
  dat <- model.frame(formula, data, na.action=na.pass)
  nn <- dim(dat)[1]
  
  if (typeof(data) == "environment") {
    id <- id
    weights <- weights
    if(is.null(call$weights)) weights <- rep(1, nn)
#    waves <- waves
  }
  else{
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
    } else if(is.null(call$weights)){
      weights <- rep.int(1,nn)
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
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves
  
  
  if(!is.numeric(dat$waves) & !is.null(dat$waves)) stop("waves must be either an integer vector or NULL")
  
  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  # Note that we need to assign weight 0 to rows with NAs
  # in order to preserve the correlation structure
  na.inds <- NULL
  
  if(any(is.na(dat))){
    na.inds <- which(is.na(dat), arr.ind=T)
  }
  
  #SORT THE DATA ACCORDING TO WAVES
  if(!is.null(waves)){
    dat <- dat[order(id, waves),]
  }else{
    dat <- dat[order(id),]
  }
  
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(corstr, cor.vec)
  
  if(is.na(cor.match)){stop("Unsupported correlation structure")}
  
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
    if( !is.element(cor.match,c(1,3)) & (sum(incomp) > 0) & !nodummy){
      dat <- dummyrows(formula, dat, incomp, maxwave, wavespl, idspl)
      id <- dat$id
      waves <- dat$waves
      weights <- dat$weights
    }
  }
  
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      if(is.factor(dat[,i])){
        dat[na.inds[,1], i] <- levels(dat[,i])[1]
      }else{
        dat[na.inds[,1], i] <- median(dat[,i], na.rm=T)
      }
    }
  }
  
  
  includedvec <- weights>0
  
  
  inclsplit <- split(includedvec, id)
  
  dropid <- NULL
  allobs <- T
  if(any(!includedvec)){
    allobs <- F
    for(i in 1:length(unique(id))){
      if(all(!inclsplit[[i]])){
        dropid <- c(dropid, unique(id)[i])
      }
    }
  }
  
  dropind <- c()
  
  if(is.element(cor.match, c(1,3))){
    dropind <- which(weights==0)
  }else if(length(dropid)>0){
    dropind <- which(is.element(id, dropid))
  }
  if(length(dropind) > 0){
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]

    id <- id[-dropind]
  }
  nn <- dim(dat)[1]
  K <- length(unique(id))
  
  
  modterms <- terms(formula)
  
  X <- model.matrix(formula,dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)
  
  p <- dim(X)[2]
  
  
  
  ### if no offset is given, then set to zero
  if(is.null(offset)){
    off <- rep(0, nn)
  }else{
    off <- offset
  }
  
  
  # Is there an intercept column?
  interceptcol <- apply(X==1, 2, all)
  
  ## Basic check to see if link and variance functions make any kind of sense
  linkOfMean <- LinkFun(mean(Y[includedvec])) - mean(off)
  
  if( any(is.infinite(linkOfMean) | is.nan(linkOfMean)) ){
    stop("Infinite or NaN in the link of the mean of responses.  Make sure link function makes sense for these data.")
  }
  if( any(is.infinite( VarFun(mean(Y))) | is.nan( VarFun(mean(Y)))) ){
    stop("Infinite or NaN in the variance of the mean of responses.  Make sure variance function makes sense for these data.")
  }
  
  if(is.null(init.beta)){
    if(any(interceptcol)){
      #if there is an intercept and no initial beta, then use link of mean of response
      init.beta <- rep(0, dim(X)[2])
      init.beta[which(interceptcol)] <- linkOfMean
    }else{
      stop("Must supply an initial beta if not using an intercept.")
    }
  }
  
  
  # Number of included observations for each cluster
  includedlen <- rep(0, K)
  len <- rep(0,K)
  uniqueid <- unique(id)
  
  tmpwgt <- as.numeric(includedvec)
  idspl <-ifelse(tmpwgt==0, NA, id)
  includedlen <- as.numeric(summary(split(Y, idspl, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
  # Get vector of cluster sizes... remember this len variable
  #len <- as.numeric(summary(split(Y, id, drop=T))[,1])
  
  # Set the initial alpha value
  if(is.null(init.alpha)){
    alpha.new <- 0.2
    if(cor.match==4){
      # If corstr = "m-dep"
      alpha.new <- 0.2^(1:Mv)
    }else if(cor.match==5){
      # If corstr = "unstructured"
      alpha.new <- rep(0.2, sum(1:(max(len)-1)))
    }else if(cor.match==7){
      # If corstr = "userdefined"
      alpha.new <- rep(0.2, max(unique(as.vector(corr.mat))))
    }
  }else{
    alpha.new <- init.alpha
  }
  #if no initial overdispersion parameter, start at 1
  if(is.null(init.phi)){
    phi <- 1
  }else{
    phi <- init.phi
  }
  
  beta <- init.beta
  
  
  #Set up matrix storage
  StdErr <- Diagonal(nn)
  dInvLinkdEta <- Diagonal(nn)
  Resid <- Diagonal(nn)
  #  if( (max(len)==1) & cor.match != 1 ){
  #    warning("Largest cluster size is 1. Changing working correlation to independence.")
  #    cor.match <- 1
  #    corstr <- "independence"
  #  }
  
  # Initialize for each correlation structure
  if(cor.match == 1){
    # INDEPENDENCE
    R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
    BlockDiag <- getBlockDiag(len)$BDiag
  }else if(cor.match == 2){
    # AR-1
    tmp <- buildAlphaInvAR(len)
    # These are the vectors needed to update the inverse correlation
    a1<- tmp$a1
    a2 <- tmp$a2
    a3 <- tmp$a3
    a4 <- tmp$a4
    # row.vec and col.vec for the big block diagonal of correlation inverses
    # both are vectors of indices that facilitate in updating R.alpha.inv
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    BlockDiag <- getBlockDiag(len)$BDiag
    
  }else if(cor.match == 3){
    # EXCHANGEABLE
    # Build a block diagonal correlation matrix for updating and sandwich calculation
    # this matrix is block diagonal with all ones.  Each block is of dimension cluster size.
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    
    #Create a vector of length number of observations with associated cluster size for each observation
    n.vec <- vector("numeric", nn)
    index <- c(cumsum(len) - len, nn)
    for(i in 1:K){
      n.vec[(index[i]+1) : index[i+1]] <-  rep(includedlen[i], len[i])
    }
    
    # n.vec <- vector("numeric", nn)
    # index <- c(cumsum(len) - len, nn)
    # for(i in 1:K){
    #  n.vec[(index[i]+1) : index[i+1]] <-  rep(len[i], len[i])
    # }
  }else if(cor.match == 4){
    # M-DEPENDENT, check that M is not too large
    if(Mv >= max(len)){
      stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
    }
    
    # Build block diagonal similar to in exchangeable case, also get row indices and column
    # indices for fast matrix updating later.
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    R.alpha.inv <- NULL
  }else if(cor.match == 5){
    # UNSTRUCTURED
    if( max(len^2 - len)/2 > length(len)){
      stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
    }
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
  }else if(cor.match == 6){
    # FIXED
    # check if matrix meets some basic conditions
    corr.mat <- checkFixedMat(corr.mat, len)
    
    R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
    BlockDiag <- getBlockDiag(len)$BDiag
  }else if(cor.match == 7){
    # USERDEFINED
    corr.mat <- checkUserMat(corr.mat, len)
    
    # get the structure of the correlation matrix in a way that
    # I can use later on.
    tmp1 <- getUserStructure(corr.mat)
    corr.list <- tmp1$corr.list
    user.row <- tmp1$row.vec
    user.col <- tmp1$col.vec
    struct.vec <- tmp1$struct.vec
    
    # the same block diagonal trick.
    tmp2 <- getBlockDiag(len)
    BlockDiag <- tmp2$BDiag
    row.vec <- tmp2$row.vec
    col.vec <- tmp2$col.vec
    
  }else if(cor.match == 0){
    stop("Ambiguous Correlation Structure Specification")
  }else{
    stop("Unsupported Correlation Structure")
  }
  
  stop <- F
  converged <- F
  count <- 0
  beta.old <- beta
  unstable <- F
  phi.old <- phi
  
  
  # Main fisher scoring loop
  while(!stop){
    count <- count+1
    
    eta <- as.vector(X %*% beta) + off
    
    mu <- InvLink(eta)
    
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    if(!scale.fix){
      phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, includedlen, sqrtW, useP)
    }
    phi.new <- phi
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(cor.match == 2){
      # AR-1
      alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len,
                                 StdErr, p, included, includedlen,
                                 includedvec, allobs, sqrtW,
                                 BlockDiag, useP)
      R.alpha.inv <- getAlphaInvAR(alpha.new, a1,a2,a3,a4, row.vec, col.vec)/phi
    }else if(cor.match == 3){
      #EXCHANGEABLE
      alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, len, StdErr,
                                 Resid, p, BlockDiag, included,
                                 includedlen, sqrtW, useP)
      #R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
      R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
    }else if(cor.match == 4){
      # M-DEPENDENT
      if(Mv==1){
        alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len,
                                   StdErr, p, included, includedlen,
                                   includedvec, allobs,
                                   sqrtW, BlockDiag, useP)
      }else{
        alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, id, len,
                                     StdErr, Resid, p, BlockDiag, Mv,
                                     included, includedlen,
                                     allobs, sqrtW, useP)
        if(sum(len>Mv) <= p){
          unstable <- T
        }
      }
      if(any(alpha.new >= 1)){
        stop <- T
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, col.vec)/phi
    }else if(cor.match == 5){
      # UNSTRUCTURED
      alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, id, len,
                                      StdErr, Resid,  p, BlockDiag,
                                      included, includedlen, allobs,
                                      sqrtW, useP)
      # This has happened to me (greater than 1 correlation estimate)
      if(any(alpha.new >= 1)){
        stop <- T
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, row.vec, col.vec)/phi
    }else if(cor.match ==6){
      # FIXED CORRELATION, DON'T NEED TO RECOMPUTE
      R.alpha.inv <- R.alpha.inv*phi.old/phi
    }else if(cor.match == 7){
      # USER SPECIFIED
      alpha.new <- updateAlphaUser(Y, mu, phi, id, len, StdErr, Resid,
                                   p, BlockDiag, user.row, user.col,
                                   corr.list, included, includedlen,
                                   allobs, sqrtW, useP)
      R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
    }else if(cor.match == 1){
      # INDEPENDENT
      R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
      alpha.new <- "independent"
    }
    
    
    beta.list <- updateBeta(Y, X, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, W, included)
    beta <- beta.list$beta
    
    phi.old <- phi
    if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){converged <- T; stop <- T}
    if(count >= maxit){stop <- T}
    beta.old <- beta
  }
  biggest <- which.max(len)[1]
  index <- sum(len[1:biggest])-len[biggest]
  
  if(K == 1){
    biggest.R.alpha.inv <- R.alpha.inv
    if(cor.match == 6) {
      biggest.R.alpha <- corr.mat*phi
    }else{
      biggest.R.alpha <- solve(R.alpha.inv)
    }
  }else{
    biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
    if(cor.match == 6){
      biggest.R.alpha <- corr.mat[1:len[biggest] , 1:len[biggest]]*phi
    }else{
      biggest.R.alpha <- solve(biggest.R.alpha.inv)
    }
  }
  
  
  eta <- as.vector(X %*% beta) + off
  if(sandwich){
    sandvar.list <- getSandwich(Y, X, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, 
                                beta.list$hess, StdErr, dInvLinkdEta, BlockDiag, W, included)
  }else{
    sandvar.list <- list()
    sandvar.list$sandvar <- "no sandwich"
  }
  
  if(!converged){warning("Did not converge")}
  if(unstable){warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")}
  
  
  # Create object of class geem with information about the fit
  dat <- model.frame(formula, data, na.action=na.pass)
  X <- model.matrix(formula, dat)
  
  if(is.character(alpha.new)){alpha.new <- 0}
  results <- list()
  results$beta <- as.vector(beta)
  results$phi <- phi
  results$alpha <- alpha.new
  if(cor.match == 6){
    results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat,1)!=0)])
  }
  results$coefnames <- colnames(X)
  results$niter <- count
  results$converged <- converged
  results$naiv.var <- solve(beta.list$hess)  ## call model-based
  results$var <- sandvar.list$sandvar
  results$call <- call
  results$corr <- cor.vec[cor.match]
  results$clusz <- len
  results$FunList <- famret
  results$X <- X
  results$offset <- off
  results$eta <- eta
  results$dropped <- dropid
  results$weights <- weights
  results$terms <- modterms
  results$y <- Y
  results$biggest.R.alpha <- biggest.R.alpha/phi
  results$formula <- formula
  class(results) <- "geem"
  return(results)
}


#####################################################################################
## Not exported below
#####################################################################################

#Do actual gee fitting
geem.fit <- function()



### Simple moment estimator of dispersion parameter
updatePhi <- function(YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW, useP){
  nn <- sum(includedlen)
  
  resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu))
  
  phi <- (1/(sum(included)- useP * p))*crossprod(resid, resid)
  
  return(as.numeric(phi))
}

### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
updateBeta = function(YY, XX, beta, off, InvLinkDeriv, InvLink,
                      VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, W, included){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(XX%*%beta.new) + off
    
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
getSandwich = function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv,
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


######git setup test delete this

