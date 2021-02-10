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
#' @param waves An integer vector identifying the time ordering within clusters 
#' (i.e. levels of \code{id}). Note that only the ordering is used, NOT the 
#' numeric values of the wave. This means that non-equidistant time points are
#' not being used for fitting in e.g. AR1. 
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
#'    
#' @param corr.mat The correlation matrix for \code{"fixed"}.  Matrix should
#'   be symmetric with dimensions >= the maximum cluster size.  If the correlation 
#'   structure is \code{"userdefined"}, then this is a matrix describing which 
#'   correlations are the same.
#'    
#' @param useP If set to \code{FALSE}, do not use the n-p correction for 
#'    dispersion and correlation estimates, as in Liang and Zeger. This can be 
#'    useful when the number of observations is small, as subtracting p may yield 
#'    correlations greater than 1.
#'  
#' @output BLABLABLA
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
#'  matrix with consecutive integers.  \code{geelm} only looks at the upper 
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
#'   as \code{(LinkFun, VarFun, InvLink, InvLinkDeriv)}, then \code{geelm} 
#'   assumes that the functions are given in that order.  LinkFun and VarFun 
#'   are the link and variance functions. InvLink and InvLinkDeriv are the inverse 
#'   of the link function and the derivative of the inverse of the link function 
#'   and so are decided by the choice of the link function.
#' 
#' @return An object of class \code{geelm} (inherits from \code{geeglm}) representing the fit.
#' 
#'   Observations with a \code{NA} in the variables specified in the \code{formula}
#'   argument, \code{weights} (if used) or \code{waves} (if used) will be assigned 
#'   a weight of 0.  Note that these weights  are now the same as PROC GEE weights 
#'   and not PROC GENMOD. CHECK IF THIS IS STILL TRUE OR SHOULD BE CHANGED. 
#' 
#' @author Anne Helby Petersen, Lee McDaniel & Nick Henderson
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
#' far1 <- geelm(count~ time, id=subject ,data = testdat, family=poisson, 
#'              corstr="ar1")
#'              
#' ### Ohio respiratory data from geepack
#'  if(require(geepack)){
#'      data("ohio", package="geepack")
#'      resplogit <- geelm(resp ~ age + smoke + age:smoke, id=id, data = ohio, 
#'                        family = binomial, corstr = "m-dep" , Mv = 1)
#'      LinkFun <- function(arg){qcauchy(arg)}
#'      InvLink <- function(arg){pcauchy(arg)}
#'      InvLinkDeriv <- function(arg){dcauchy(arg)}
#'      VarFun <- function(arg){arg*(1-arg)}
#'      FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)
#'      
#'      respcauchit <- geelm(resp ~ age + smoke + age:smoke, id=id, data = ohio, 
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
#'     seiz <- geelm(y~ x + trt + x:trt+ offset(log(t)), id=id,data = seiz.l, 
#'                  family = poisson, corstr = "exchangeable")
#' }
#' 
#' @export
geelm <- function(formula, id, waves=NULL, data = parent.frame(),
                 family = gaussian, corstr = "independence", Mv = 1,
                 weights = NULL, corr.mat = NULL, 
                 offset = NULL,
                 nodummy = FALSE,  
                 output = "geelm",
                 control = geelm.control()){
  
  ########################################################################
  #Check and prep input arguments ########################################
  ########################################################################
  
  thiscall <- match.call()
  
  
  
  ### First, get all the relevant elements from the arguments
  dat <- model.frame(formula, data, na.action = na.pass)
  nn <- dim(dat)[1]
  
  #Make id / weight / waves argument so that they can match character
  #strings OR names provided in enclosing environment/data
  #??? doesn't work with character strings??? 
  if(typeof(data) != "environment") {
    if(length(thiscall$id) == 1){
      subj.col <- which(colnames(data) == thiscall$id)
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      } else {
        id <- eval(thiscall$id, envir=parent.frame())
      }
    } else if(is.null(thiscall$id)) {
      id <- 1:nn
    }
    
    if(length(thiscall$weights) == 1) {
      weights.col <- which(colnames(data) == thiscall$weights)
      if(length(weights.col) > 0) {
        weights <- data[,weights.col]
      } else {
        weights <- eval(thiscall$weights, envir=parent.frame())
      }
    } 
    
    if(length(thiscall$waves) == 1) {
      waves.col <- which(colnames(data) == thiscall$waves)
      if(length(waves.col) > 0) {
        waves <- data[,waves.col]
      } else {
        waves <- eval(thiscall$waves, envir=parent.frame())
      }
    } else if(is.null(thiscall$waves)) {
      waves <- NULL
    }
  }
  
  # Initialize weights if not supplied by user
  if (is.null(weights)) weights <- rep(1, nn)
  
  # Check waves argument
  if (!is.null(waves) && !identical(round(waves, 0), waves)) stop("waves must be either an integer vector or NULL") 
  
  
  # Store objects for output
  prior.weights <- weights
  
  
  # Organize all dataset in dat
  dat$id <- id
  dat$weights <- weights
  dat$waves <- waves #!!! this may be NULL => not controlled what slots are available in dat

  
  # Sort data. If waves are available, sort accord to id and waves, otherwise 
  # sort only according to id
  if (!is.null(waves)) {
    neworder <- order(id, waves)
  } else {
    neworder <- order(id)
  }
  dat <- dat[neworder, ] 
  
  # Check if waves are equidistant. If they are, we are done using the 
  # values of waves now, as data have already been sorted according to waves within ids.
  # If the waves are not equidistant, signal a warning and proceed as if they were 
  # equidistant (i.e. use only order information - this has already been done in sorting)
  if (!is.null(waves)) {
    waves_are_equidist <- by(dat$waves, as.factor(dat$id), is_equidistant)
    if (any(!waves_are_equidist)) {
      warning(paste("Non-equidistant waves were provided.",
                    "Note that only their ordering was used for model fitting.",
                    "Their numeric values were ignored."))
    }
  }
  
  
  # Find missing information in variables used for the linear predictor and response
  # note: na.inds contains information on both row and columns of NAs - both are needed below. 
  na.inds <- NULL
  if(any(is.na(dat))){
    na.inds <- which(is.na(dat), arr.ind = TRUE)
  }
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", 
               "fixed", "userdefined")
  corstr <- cor.vec[charmatch(corstr, cor.vec)]
  
  #!!!to do:
  # - add step to change corstr if settings simplify (for example "m-dependent" and Mv == 1 => "ar1"?)
  
  #!!! to do:
  #!!! check what happens if length(cor.str) > 1 and check what happens if there is no match (old code below
  # NOT replaced yet):
  #*#  if(is.na(cor.match)){stop("Unsupported correlation structure")} #!!check if this still works
  #*#
  #*#   else if(cor.match == 0){
  #*#         stop("Ambiguous Correlation Structure Specification")
  #*#   }else{ 
  #*#   stop("Unsupported Correlation Structure")
  #*#   }
  
  

  # Note that we need to assign weight 0 to rows with NAs
  # in order to preserve the correlation structure
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0 #!!! consider: should this be done in "weights" within dat instead?
  } 
  
  includedvec <- weights > 0
  inclsplit <- split(includedvec, id)
  
  #Find indexes to be dropped
  dropind <- NULL 
  
  # Identify ids that are missing for all observations, i.e. whole cluster must be 
  # missing. If so, omit it, otherwise ignore. 
  if(corstr %in% c("independence", "exchangeable")) {
    dropind <- which(weights == 0)
  } else {
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
        }
      }
      if (sum(dodropid) > 0) {
          dropind <- which(id %in% uniqueid[dodropid])  
      }
    }
  }
  #drop indexes just found
  if(length(dropind) > 0) {
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- id[-dropind]
  }

    
  X <- model.matrix(formula, dat) #nas 
  Y <- model.response(dat)
  
  #Handle offset. Note: Can be specified both in formula AND
  #in seperate offset argument (needed for glm type methods including
  #anova comparisons of nested models)
  ## if no offset is given, then set to zero
  formulaoffset <- model.offset(dat)
  if (is.null(formulaoffset)) formulaoffset <- rep(0, nrow(X))
  
  if (is.null(offset)) {
    offset <- formulaoffset
  } else {
    offset <- offset + formulaoffset
  }
  
 
  
  
  # add extra info to corstr if necessary. These are stored as attributes. 
  if (corstr == "m-dependent") {
    attr(corstr, "Mv") <- Mv
  } else if (corstr %in% c("fixed", "userdefined")) {
    attr(corstr, "corr.mat") <- corr.mat
  }
  
  # handle family argument
  famret <- getfam(family)
  
  #rename function names for standardized output family object.
  # note: not sure why original geeM had different function names
  # for user-supplied families. May be possible to drop this
  # convention altogether
  #browser()
  if (!inherits(famret, "family")) {
    names(famret)[c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv")] <-
      c("linkfun", "variance", "linkinv", "mu.eta")
  }
  

  ##########################################################################################
  # Do actual fitting      #################################################################
  ##########################################################################################
  
  #!!! check: can we avoid passing allobs argument and instead recompute it 
  # in geem.fit?
  
  results <- geelm.fit(x = X, y = Y, offset = offset, weights = weights,
                  control = control, id = id, family = famret,
                  corstr = corstr, allobs = allobs)
  

  
  ##########################################################################################
  # Pack and return output #################################################################
  ##########################################################################################
  
  if (output == "geeM") {
    # Create object resembling original geeM::geem() output

    #add slots not already available in geelm.fit output
    results$coefnames <- colnames(X)
    results$call <- thiscall
    results$X <- X
    results$dropped <- uniqueid[dodropid]
    results$terms <- terms(formula)
    results$y <- Y
    results$formula <- formula
    results$var <- results$vbeta
    
    # new dat and X objects are standard in geem - these will be ordered as 
    # the original input data but may differ in observations if NAs were dropped
    # along the way!
    newdat <- model.frame(formula, data, na.action = na.pass) 
    results$X <- model.matrix(formula, newdat)

    #reorder list to make it identical to previous structure
    old_geem_out_order <- c("beta", "phi", "alpha", "coefnames", 
                            "niter", "converged", "naiv.var", "var",
                            "call", "corr", "clusz", "FunList", 
                            "X", "offset", "eta", "weights", "terms",
                            "y", "biggest.R.alpha", "formula")
    results <-results[old_geem_out_order]
    
    class(results) <- "geem"
    return(results)
  } 
  
  if (output %in% c("geeglm", "geeglm")) {
#    browser()
    coefs <- results$beta
    names(coefs) <- colnames(X)
    
    # Rank of model matrix
    model_rank <- Matrix::rankMatrix(na.omit(X)) #!!!!! check: ok to do na.omit here?
    
    # construct modelinfo which is used both for geese and full geeglm object
    modelInfo = list(mean.link  = famret$link,
                     variance = famret$family,
                     sca.link = "identity",
                     cor.link = "identity",
                     corstr = corstr,
                     scale.fix = control$scale.fix)
   
    # construct variance objects with correct dimensions but filled with zeros
    # as these alternative variance estimation options are not yet supported
    # note: zero matrices in this scenario is geepack standard 
    vbeta <- results$vbeta #note: need "regular" matrix for geeglm methods to work
    vbeta_otherse <- vbeta
    vbeta_otherse[,] <- 0
    vgamma_otherse <- matrix(0, nrow = length(results$phi), ncol = length(results$phi))
    valpha_otherse <- matrix(0, nrow = length(results$alpha), ncol = length(results$alpha))
    
    # construct geese object - necessary for geepack methods
    geeseobj <- list(beta = coefs,
                     vbeta = vbeta,
                     vbeta.j1s = vbeta_otherse,
                     vbeta.fij = vbeta_otherse,
                     vbeta.ajs = vbeta_otherse,
                     vbeta.naiv = as.matrix(results$naiv.var),
                     gamma = results$phi, #!!! OBS: Tjek om phi og gamma er det samme eller om der er en transformation imellem
                     vgamma =  vgamma_otherse, #!! placeholder - check: do we compute this at all?
                     vgamma.j1s = vgamma_otherse,
                     vgamma.fij = vgamma_otherse,
                     vgamma.ajs = vgamma_otherse,
                     alpha = results$alpha, #!! placeholde - check: do we compute this at all?
                     valpha =  valpha_otherse,
                     valpha.j1s = valpha_otherse,
                     valpha.fij = valpha_otherse,
                     valpha.ajs = valpha_otherse,
                     model = modelInfo,
                     control = control,
                     error = NA,  #not sure what this one does
                     clusz = results$clusz,
                     zsca.names = NULL, #add printed name for dispersion parameter here
                     zcor.names = NULL, #add printed names for alpha parameters here
                     xnames = colnames(X)
                     )
    class(geeseobj) <- c("geese", "list")
    
    #leave out dropped observations from ordering 
    if (length(dropind) > 0) { #case: any obs dropped 
      oldorder_noNA <- order(neworder[-dropind])
    } else { #case: no obs dropped
      oldorder_noNA <- order(neworder)
    }
  
    out <- list(coefficients = coefs,
                    residuals = results$resid[oldorder_noNA],
                    fitted.values = results$fitted.values[oldorder_noNA],
                    effects = NA, #not required for glm objects, ever used? note: documented in white book
                    rank = model_rank, #rank of model matrix
                    qr = qr(na.omit(X)), #not entirely sure if this is the correct matrix to compute QR decompostion for... 
                                         #also: should we reinset dropped rows and fill them with NAs post hoc?
                                         #also: order back to oldorder?
                    family = famret,
                    linear.predictors = results$eta[oldorder_noNA],
                    weights = results$weights[oldorder_noNA],
                    prior.weights = prior.weights,
                    df.residual = sum(results$weights != 0) - model_rank,
                    y = Y[oldorder_noNA],
                    model = dat[oldorder_noNA,], #note: this is reordered (back to input order) in order to get correct
                                                 #outout from anova.
                    call = thiscall,
                    formula = formula,
                    terms =  terms(formula),
                    data = data,
                    offset = results$offset[oldorder_noNA],
                    control = control,
                    method = "geem.fit",
                    constrasts = attr(X, "contrasts"),
                    xlevels = get_xlevels(dat),
                    geese = geeseobj, 
                    modelInfo = modelInfo,
                    id = id,
                    corstr = corstr,
                    cor.link = "identity",
                    std.err = control$std.err)
    class(out) <- c("geelm", "geeglm", "gee", "glm", "lm")
    
    return(out)
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
  conv <- FALSE
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



# Get levels for all factors from model.frame type object
get_xlevels <- function(modelframe) {
  classes <- attr(terms(modelframe), "dataClasses")
  factors <- names(classes[classes == "factor"])
  nfactors <- length(factors)
  if (nfactors > 0) { 
    xlevels <- as.list(1:nfactors) 
    for (i in 1:nfactors) {
      xlevels[[i]] <- levels(modelframe[, factors[i]])
    }
    names(xlevels) <- factors
    return(xlevels)
  }
  NULL
}


# Check if a sorted numeric vector, x, is equidistant, i.e. same difference
# for all two subsequent entries. If x contains one or no elements, TRUE
# is outputted
is_equidistant <- function(x) {
  nx <- length(x)
  if (nx > 1) {
    return(length(unique(x[-1] - x[-nx])) == 1)
  } else {
    return(TRUE)
  }
}

