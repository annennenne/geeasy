#' Fit Generalized Estimating Equation-based Linear Models
#' 
#' Estimate mean structure parameters and their corresponding standard errors for
#' generalized linear models with clustered or correlated observations by use of 
#' generalized estimating equations.  
#' 
#' @inheritParams stats::glm
#' 
#' @param formula A formula expression similar to that for \code{\link{glm}}, 
#  of the form \code{response~predictors}.  An offset is allowed, as in \code{glm}.
#' 
#' @param id A vector identifying the clusters. If NULL, then each observation is 
#' assigned its own cluster.
#' 
#' @param waves An numeric vector identifying the time ordering within clusters 
#' (i.e. levels of \code{id}). By default, data are assumed 
#' to be sorted such that observations in a cluster are in consecutive rows 
#' and higher numbered rows in a cluster are assumed to be later. Note that only the 
#' ordering of the values in \code{waves} is used, NOT the numeric values themselves. 
#' This means that e.g. having waves equal to \code{c(1, 2, 3)} 
#' or \code{c(1, 2, 7)} within a cluster results in the same model.
#' 
#' @param data An optional data frame containing the variables in the model.
#' 
#' @param family A description of the error distribution and link function to be used 
#' in the model. The argument can be one of three options: a \code{family} object, 
#' a character string, or a list of functions. For more information on how to use \code{family} 
#'  objects, see Details below. 
#'  
#' @param corstr A character string specifying the correlation structure. 
#' The default is "independence". Allowed structures are: \code{"independence"}, 
#' \code{"exchangeable"},  \code{"ar1"}, \code{"m-dependent"}, \code{"unstructured"}, 
#' \code{"fixed"}, and \code{"userdefined"}.  Any unique substring may be supplied.  
#'  If \code{"fixed"} or \code{"userdefined"}, then \code{corr.mat} must be 
#'  specified.  If \code{"m-dependent"}, then \code{Mv} is relevant.
#'  
#' @param Mv For \code{"m-dependent"}, the value for \code{m}. 
#'    
#' @param corr.mat The correlation matrix for \code{"fixed"}.  Matrix should
#'   be symmetric with dimensions >= the maximum cluster size.  If the correlation 
#'   structure is \code{"userdefined"}, then this is a matrix describing which 
#'   correlations are the same.
#'   
#' @param control A list of parameters for controlling the fitting process. 
#'    
#' @param output Output object type. There are two options; 1) \code{"geelm"} (default), resulting in 
#' an output that inherits the structure of \code{geepack}s \code{geeglm} object, or 2)
#' \code{"geem"} (or its alias \code{"geeM"}) which results in an output that has the structure
#'  of \code{geeM}s \code{geem} object. 
#' 
#' @param engine Engine used to fit the model. The default, \code{"geeasy"} uses this
#' package (built on the \code{geeM} package), while \code{"geepack"} uses
#' the function \code{geeglm} from \code{geepack} to fit the model. Note that if 
#' the geepack engine is used, the data are sorted according to id (and possibly 
#' waves within id) and NAs are dropped before the data is used 
#' (this differs from the standard in geepack).
#'     
#' @details Users may specify functions for link and variance functions, but the
#'  functions must be vectorized functions. 
#'    
#'  Offsets can be specified in the model formula, as in \code{glm()} or they may be
#'  specified using the \code{offset} argument. If offsets are specified in both ways, 
#'  their sum is used as an offset.
#'  
#'  For the \code{"userdefined"} correlation option, the function accepts a 
#'  matrix with consecutive integers. Each such integer represent a distinct
#'  parameter that will be estimated.  All entries given as 1 will be assumed
#'   to be the same as each other and will be assumed to be possibly different 
#'   from entries with a 2, and so on.\code{geelm} only looks at the upper 
#'  triangle of the matrix.  Any entry given as 0 will be fixed at 0.  
#'   
#'  If observations are dropped because they have a weight of 0, then the 
#'  denominator for the moment estimates of the correlation matrices are 
#'  calculated using the number of non-zero Pearson residuals for the 
#'  correlation structures \code{unstructured}, \code{userdefined} and 
#'  \code{m-dependent} with \code{Mv>1}.  Therefore, residuals numerically 
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
#' It contains the following slots:
#' 
#' \code{$coefficients}: Coefficients from the mean structure model (betas) on their 
#' original scales
#' 
#' \code{$residuals}: Pearson residuals, in the order of the inputted dataset (with NAs omitted).
#' 
#' \code{$fitted.values}: Fitted values (response scale), in the order of the inputted dataset 
#'  (with NAs omitted).
#'  
#' \code{$rank}: The rank of the model matrix, i.e. the number of estimated mean structure
#' coefficients.
#' 
#' \code{$qr}: QR decomposition of the model matrix (NA omitted). 
#' 
#' \code{$family}: A family object specifying which exponential family was used for fitting
#' the mean structure model, see \code{\link{family}} for more information. 
#' 
#' \code{$linear.predictors}: The linear predictor on the original scale.
#' 
#' \code{$weights}: Weights used for computations, in the order of the inputted dataset
#'  (NAs omitted).
#'  
#' \code{$prior.weights}: The original weights used to produce this geeglm object (set
#' by user or defaulted to 1 for all observations).
#' 
#' \code{$df.residuals}: Residual degrees of freedom.
#' 
#' \code{$y}: Outcome variable, in the order of the inputted dataset (NAs omitted).
#' 
#' \code{$model}: The model.frame, ordered as the original inputted data with NAs omitted.
#' 
#' \code{$call}: The original function call that produced this geeglm object.
#' 
#' \code{$formula}: The formula used in the original call.
#' 
#' \code{$terms}: The terms of the formula used in the original call.
#' 
#' \code{$data}: The original dataset that was used for producing this geeglm object.
#' 
#' \code{$offset}: Offset used for fitting the model, ordered as the original inputted data
#' with NAs omitted.
#' 
#' \code{$control}: Value of control parameters used for fitting the model. 
#' 
#' \code{$method}: Internal function used for fitting the model.
#' 
#' \code{$contrasts}: Contrasts used in the model matrix.
#' 
#' \code{$xlevels}: Levels of factor variables used in the model formula (if any).
#' 
#' \code{$geese}: An object containing further information about the variance estimation, 
#' including a variance matrix for the beta-coefficients (\code{$vbeta}), the estimated 
#' coefficients for the working correlation matrix (\code{$alpha}), the estimated dispersion 
#' parameter (\code{$gamma}), and the individual cluster sizes (\code{$clusz}). See 
#' \code{\link{geese}} for more information. 
#' 
#' \code{$modelInfo}: Information about the link functions used for fitting the mean, variance 
#' and scale structures of the model. 
#' 
#' \code{$id}: IDs used for identifying the clusters, ordered as the original inputted data 
#' with NAs omitted.
#' 
#' \code{$corstr}: Name of the correlation structured imposed on the model. If the 
#' correlation structure requires further information, it is stored in a suitably named
#' attribute. For example, for m-dependent correlation structures, the m scalar is available
#' in an attribute named \code{Mv}. 
#' 
#' \code{$cor.link}: Link function used for the correlation structure.
#' 
#' \code{$std.err}: Method used to estimate the standard error of the mean structure 
#' coefficients (betas).
#' 
#' @author Anne Helby Petersen, Lee McDaniel & Nick Henderson
#' 
#' @seealso \code{\link{glm}}, \code{\link{formula}}, \code{\link{family}}
#' 
#' @keywords models robust
#' 
#' @importFrom stats complete.cases gaussian model.frame model.offset model.response na.omit na.pass 
#' @examples
#' 
#' # load data
#' data("respiratory")
#' respiratory$useid <- interaction(respiratory$center, respiratory$id)
#' 
#' # fit model
#' m <- geelm(outcome ~ treat + sex + age + baseline, 
#'            data = respiratory, id = useid,
#'                       family = "binomial", corstr = "exchangeable")
#' 
#' \dontrun{
#' get_jack_se <- function(object, dat){
#'     parm <- sapply(1:nrow(dat),
#'                    function(i){
#'                        dat.i <- dat[-i,]
#'                        coef(update(object, data=dat.i))
#'                    })
#'     parm <- t(parm)
#'     parm.mean <- apply(parm, 2, mean)
#'     
#'     parm.cent <- sapply(1:nrow(parm),
#'                         function(i){
#'                             parm[i, ] - parm.mean
#'                         })
#'     parm.cent <- t(parm.cent) 
#'     
#'     jack.var <- ((nrow(dat)-1) / nrow(dat)) * t(parm.cent) %*% parm.cent
#'     jack.se <- sqrt(diag(jack.var))
#'     jack.se
#' }
#' 
#' 
#' # load data
#' data("respiratory")
#' respiratory$useid <- interaction(respiratory$center, respiratory$id)
#' 
#' # fit model
#' obj <- geelm(outcome ~ treat + sex + age + baseline, 
#'            data = respiratory, id = useid,
#'                       family = "binomial", corstr = "exchangeable")
#' 
#' dat <- respiratory
#' get_jack_se(obj, dat)
#' summary(obj) |> coef()
#' }
#' 
#' @importFrom geepack geeglm geese.control
#' 
#' @export
geelm <- function(formula, id = NULL, waves = NULL, data = parent.frame(),
                 family = gaussian, corstr = "independence", Mv = 1,
                 weights = NULL, corr.mat = NULL, 
                 offset = NULL,
                 engine = "geeasy",
                 output = "geelm",
                 control = geelm.control()
                 ){
  
  #########################################################################################
  #Check and handle input arguments #######################################################
  #########################################################################################
  
  thiscall <- match.call()
  
  ########################################################################
  # Handle correlation structure arguments
  ########################################################################
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", 
               "unstructured", "fixed", "userdefined")
  corstr <- cor.vec[charmatch(corstr, cor.vec)]
  if (length(corstr) == 0) {
    stop("Ambiguous correlation structure specification")
  } else if (is.na(corstr)) {
    stop("Unsupported correlation structure")
  } 
  
  # Add extra info to corstr if necessary. These are stored as attributes. 
  if (corstr == "m-dependent") {
    if (!is.numeric(Mv) || !(length(Mv) == 1) || round(Mv, 0) != Mv || Mv < 0) {
      stop("Mv must be an integer")
    }
    attr(corstr, "Mv") <- Mv
  } else if (corstr %in% c("fixed", "userdefined")) {
    if (!is.matrix(corr.mat)) {
      stop(paste("A matrix must be supplied as corr.mat for",
                 corstr, "correlation structure."))
      #note: further validity checking is done by geelm.fit()
    }
    attr(corstr, "corr.mat") <- corr.mat
  }
  
  ########################################################################
  # Handle family arguments
  ########################################################################

  famret <- getfam(family)
  
  #rename function names for standardized output family object.
  # note: not sure why original geeM had different function names
  # for user-supplied families. May be possible to drop this
  # convention altogether
  if (!inherits(famret, "family")) {
    names(famret)[c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv")] <-
      c("linkfun", "variance", "linkinv", "mu.eta")
  }
  
  
  ########################################################################
  # Handle data & variable arguments
  ########################################################################
  
  dat <- model.frame(formula, data, na.action = na.pass)
  nn <- dim(dat)[1]
  
  #Make id/weight/waves/offset argument so that they match in the following
  #prioritized order:
  #1) Look in data (as name)
  #2) Look in calling environment
  #3) Look in data (as character string)
  
    # Look up id. If none are supplied, assign unique id to each 
    # observation
    if (!is.null(thiscall$id)) {
      dat$id <- lookup_var_arg(thiscall$id, data, "id")
    } else {
      dat$id <- 1:nn
    }
  
    # Look up weights. If none are supplied, assign 1 to each 
    #observation
    if (!is.null(thiscall$weights)) {
      dat$weights <- lookup_var_arg(thiscall$weights, data, "weights")
      if (!is.numeric(dat$weights)) stop("weights must be numeric or NULL") 
    } else {
      dat$weights <- rep(1, nn)
    }
  
    # Look up waves. If waves are supplied, check that they are numeric. 
    # If no waves are supplied, don't do anything. 
    # Note: if waves are NULL nothing is stored in data. This is 
    # unproblematic as waves slot in data is only accessed when 
    # use_waves is TRUE. 
    if (!is.null(thiscall$waves)) {
      dat$waves <- lookup_var_arg(thiscall$waves, data, "waves")
      if (!is.numeric(dat$waves)) stop("waves must be numeric or NULL") 
      use_waves <- TRUE
    } else {
      use_waves <- FALSE
    }
  
  
    #Handle offset. Note: Can be specified both in formula AND
    #in seperate offset argument (needed for glm type methods including
    #anova comparisons of nested models)
      
      #First, formula offset: 
      formulaoffset <- model.offset(dat)
      if (is.null(formulaoffset)) formulaoffset <- rep(0, nn)
      
      #Then look up  offset from argument:
      if (!is.null(thiscall$offset)) {
        argoffset <- lookup_var_arg(thiscall$offset, data, "offset")
      } else {
        argoffset <- rep(0, nn)
      }
      
      #Finally, add the two
      dat$offset <- formulaoffset + argoffset
 
  # Store objects that we will need for final output
  prior.weights <- dat$weights
  
  
  #########################################################################################
  # Prep data #############################################################################
  #########################################################################################
  
  ########################################################################
  # Sort data & apply waves
  ########################################################################
  
  # Sort data. If waves are available, sort accord to id and waves, otherwise 
  # sort only according to id
  if (use_waves) {
    neworder <- order(dat$id, dat$waves)
  } else {
    neworder <- order(dat$id)
  }
  dat <- dat[neworder, ] 
  
  # Check if waves are equidistant. If they are, we are done using the 
  # values of waves now, as data have already been sorted according to waves within ids.
  # If the waves are not equidistant, signal a warning and proceed as if they were 
  # equidistant (i.e. use only order information - this has already been done in sorting)
  if (use_waves) {
    waves_are_equidist <- by(dat$waves, as.factor(dat$id), is_equidistant)
    if (any(!waves_are_equidist)) {
      warning(paste("Non-equidistant waves were provided.",
                    "Note that only their ordering was used for model fitting.",
                    "Their numeric values were ignored."))
    }
    
    #remove waves from data: they have served their purpose, no need to carry them
    #around any longer 
    dat$waves <- NULL
  }
  
  
  ########################################################################
  # Handle missing information
  ########################################################################
  
  # store info about dropped observations for output reordering etc
  dropind <- which(!complete.cases(dat))
  dropid <- dat$id[dropind]
  
  # drop NAs
  dat <- na.omit(dat)

  ########################################################################
  # Make model matrix and response
  ########################################################################
  
  X <- model.matrix(formula, dat) 
  Y <- model.response(dat)
  
  
  #########################################################################################
  # Do actual fitting      ################################################################
  #########################################################################################
 
  ########################################################################
  # Engine: "geepack"
  ########################################################################
  
  # call geepack::geeglm and return its output
  if (engine == "geepack") {
    
    # order data & id and drop NAs if any
    ord_data <- data[neworder, ]
    ord_argoffset <- argoffset[neworder]
    if (length(dropind) > 0) {
      ord_data <- ord_data[-dropind, ]
      ord_argoffset <- ord_argoffset[-dropind]
    }

    geepack_control <- geese.control(maxit = control$maxit,
                                              epsilon = control$tol)
    
    results <- do.call(geeglm, list(formula = formula, family = famret,
                              data = ord_data, 
                              weights = dat$weights,
                              offset = ord_argoffset,
                              control = geepack_control,
                              id = dat$id,
                              waves = dat$waves,
                              corstr = corstr,
                              scale.fix = control$scale.fix))
    results$call <- thiscall
    return(results)
  }
  
  ########################################################################
  # Engine: "geeasy"
  ########################################################################

  results <- geelm.fit(x = X, y = Y, offset = dat$offset, weights = dat$weights,
                  control = control, id = dat$id, family = famret,
                  corstr = corstr)

  
  ##########################################################################################
  # Pack and return output #################################################################
  ##########################################################################################
  
  ########################################################################
  # Output: geeM
  ########################################################################
  
  # Create object resembling original geeM::geem() output
  if (output %in% c("geem", "geeM")) {

        #add slots not already available in geelm.fit output
    results$coefnames <- colnames(X)
    results$call <- thiscall
    results$X <- X
    results$dropped <- dropid
    results$terms <- terms(formula)
    results$y <- Y
    results$formula <- formula
    results$var <- results$vbeta
    
    # new dat and X objects are standard in geem - these will be ordered as 
    # the original input data but may differ in observations if NAs were dropped
    # along the way!
    newdat <- model.frame(formula, data, na.action = na.pass) 
    results$X <- model.matrix(formula, newdat)

    #reorder list to make it identical to geem structure
    old_geem_out_order <- c("beta", "phi", "alpha", "coefnames", 
                            "niter", "converged", "naiv.var", "var",
                            "call", "corr", "clusz", "FunList", 
                            "X", "offset", "eta", "dropped", "weights", "terms",
                            "y", "biggest.R.alpha", "formula")
    results <-results[old_geem_out_order]
    
    class(results) <- "geem"
    return(results)
  } 
  
  ########################################################################
  # Output: geelm (note: mathces geeglm as well but this is undocumented)
  ########################################################################
  
  if (output %in% c("geelm", "geeglm")) {
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
                     vbeta.j1s = vbeta_otherse, #placeholder
                     vbeta.fij = vbeta_otherse, #placeholder
                     vbeta.ajs = vbeta_otherse, #placeholder
                     vbeta.naiv = as.matrix(results$naiv.var),
                     gamma = results$phi, 
                     vgamma =  vgamma_otherse, #placeholder
                     vgamma.j1s = vgamma_otherse, #placeholder
                     vgamma.fij = vgamma_otherse, #placeholder
                     vgamma.ajs = vgamma_otherse, #placeholder
                     alpha = results$alpha, #placeholder
                     valpha =  valpha_otherse, #placeholder
                     valpha.j1s = valpha_otherse, #placeholder
                     valpha.fij = valpha_otherse, #placeholder
                     valpha.ajs = valpha_otherse, #placeholder
                     model = modelInfo,
                     control = control,
                     error = NA,  #not sure what this one does
                     clusz = results$clusz,
                     zsca.names = NULL, #add printed name for dispersion parameter here
                     zcor.names = NULL, #add printed names for alpha parameters here
                     xnames = colnames(X)
                     )
    class(geeseobj) <- c("geese", "list")
    
    # leave out dropped observations from ordering 
    if (length(dropind) > 0) { #case: any obs dropped 
      oldorder_noNA <- order(neworder[-dropind])
    } else { #case: no obs dropped
      oldorder_noNA <- order(neworder)
    }
    
  
    out <- list(coefficients = coefs,
                    residuals = results$resid[oldorder_noNA],
                    fitted.values = results$fitted.values[oldorder_noNA],
                    rank = model_rank, #rank of model matrix
                    qr = qr(X), 
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
                    method = "geelm.fit",
                    constrasts = attr(X, "contrasts"),
                    xlevels = get_xlevels(dat),
                    geese = geeseobj, 
                    modelInfo = modelInfo,
                    id = dat$id[oldorder_noNA],
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

#Find variable names or expressions using variable names or plain vectors
#in data or calling environment
lookup_var_arg <- function(arg, data, name = NULL) {
  #First: Look for arg in data, and if not found there,
  #look in calling environment (grandparent frame of this function, hence 
  # n = 2)
  out <- eval(arg, envir = data, 
              enclos = parent.frame(n = 2)) 
  
  #If arg were not found in data/calling env, out will now 
  #just be a character string with the contents of arg. Hence it will 
  #have length 1. In this case, we will try to look for arg among the 
  #variable names in the data
  if (length(out) == 1) {
    out <- get(arg, envir = as.environment(data))
  }
  
  #If we still haven't found a variable in the data, raise error
  if (length(out) != nrow(data)) {
    stop(paste("Invalid"), ifelse(is.null(name), "argument", name))
  }
  
  out
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



### Simple moment estimator of dispersion parameter
#' @import Matrix
update_phi <- function(YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW, useP){
  nn <- sum(includedlen)
  
  resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu))
  
  phi <- (1/(sum(included)- useP * p))*crossprod(resid, resid)
  
  return(as.numeric(phi))
}

  
### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
update_beta <- function(YY, XX, beta, offset, InvLinkDeriv, InvLink,
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
    
    update <- solve(hess, esteq)
    
    
    beta.new <- beta.new + as.vector(update)
    
  }
  return(list(beta = beta.new, hess = hess))
}


### Calculate the sandiwch estimator as usual.
get_sandwich <- function(YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv,
                        InvLink, VarFun, hessMat, StdErr, dInvLinkdEta,
                        BlockDiag, W, included){
  
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
  scoreDiag <- Diagonal(x= YY - mu)
  BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
  
  numsand <- as.matrix(crossprod( StdErr %*% dInvLinkdEta %*% XX,  
                                  R.alpha.inv %*% W %*% StdErr %*% 
                                    BlockDiag %*% StdErr %*% W %*% 
                                    R.alpha.inv %*%  StdErr %*% 
                                    dInvLinkdEta %*% XX))
  
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  
  return(list(sandvar = sandvar, numsand = numsand))
}

