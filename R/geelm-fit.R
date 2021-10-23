#' @describeIn geelm
#'
#' @importFrom stats glm.fit
#' @export
geelm.fit <- function(x, y, id, offset, family, weights, control, corstr,
                      start = NULL) {
  #note: "start" argument included for the sake of anova.geeglm(). It is not
  #used in any way

  # Unpack family functions
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
  
  # Number of covariates (p) and observations in total
  p <- ncol(x)
  nn <- nrow(x) 
  
  ## Basic check to see if link and variance functions make any kind of sense
  meanY <- mean(y)
  linkOfMean <- LinkFun(meanY) - mean(offset)
  if(!is.finite(linkOfMean)) {
    stop("Infinite or NaN in the link of the mean of responses. Make sure link function makes sense for these data.")
  }
  if(!is.finite(VarFun(meanY))) {
    stop("Infinite or NaN in the variance of the mean of responses. Make sure variance function makes sense for these data.")
  }
  
  # Check if we are using sandwich estimation
  if (control$std.err == "san.se") {
    sandwich <- TRUE
  } else {
    sandwich <- FALSE
  }
  
  #This argument was used to handle partially missing obsevations in geem and is
  #needed for alpha-update functions. If NA handling is rewritten, consider if 
  #it should be used for anything. 
  allobs <- TRUE
  
  ###############################################################################################
  # Initialization #############################################################################
  ###############################################################################################
  
  ####################################################################
  # Initialize estimation vars #######################################
  ####################################################################
  
  #Initialize beta
  #First choice: Values supplied via control
  #Second choice: Values from glm fit
  beta <- control$init.beta
  if (is.null(beta)) {
    m_glm <- glm.fit(x = x, y = y, offset = offset, family = family, 
                     weights = weights)
    beta <- m_glm$coefficients
  }
  
  #Initialize phi
  phi <- control$init.phi
  
  #Maximum number of iterations
  maxit <- control$maxit
  
  #Additional control arguments
  scale.fix <- control$scale.fix
  useP <- control$useP  
  tol <- control$tol
  
  ####################################################################
  # Perform computations needed no matter which corstr is used #######
  ####################################################################
  
  # Number of included observations for each cluster
  uniqueid <- unique(id)
  K <- length(uniqueid)
  
  
  # includedvec is logical vector with T if weight > 0, F otherwise
  includedvec <- weights > 0 
  
  #number of complete (non-NA) obs in each cluster
  includedlen <- as.numeric(summary(split(y, ifelse(!includedvec, NA, id), 
                                          drop = TRUE))[,1])  
  
  #number of total obs in each cluster
  len <- as.numeric(summary(split(y, id, drop = TRUE))[,1]) 
  
  # W is a diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  W <- Diagonal(x = weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x = as.numeric(includedvec)) 
  
  #Set up matrix storage
  StdErr <- Diagonal(nn)
  dInvLinkdEta <- Diagonal(nn)
  Resid <- Diagonal(nn) 
  
  
  ###############################################################################################
  # Initialization according to corstr ##########################################################
  ###############################################################################################
  
  # Initialize for each correlation structure
  if(corstr == "independence") {
    R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
    BlockDiag <- get_block_diag(len)$BDiag
  } else if(corstr == "ar1"){
    tmp <- build_alpha_inv_ar(len)
    # These are the vectors needed to update the inverse correlation
    a1 <- tmp$a1
    a2 <- tmp$a2
    a3 <- tmp$a3
    a4 <- tmp$a4
    # row.vec and col.vec for the big block diagonal of correlation inverses
    # both are vectors of indices that facilitate in updating R.alpha.inv
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    BlockDiag <- get_block_diag(len)$BDiag
    
  } else if(corstr == "exchangeable"){
    # Build a block diagonal correlation matrix for updating and sandwich calculation
    # this matrix is block diagonal with all ones.  Each block is of dimension cluster size.
    tmp <- get_block_diag(len)
    BlockDiag <- tmp$BDiag
    
    #Create a vector of length number of observations with associated cluster size for each observation
    n.vec <- vector("numeric", nn)
    index <- c(cumsum(len) - len, nn)
    for(i in 1:K){
      n.vec[(index[i]+1) : index[i+1]] <-  rep(includedlen[i], len[i])
    }
    
  } else if(corstr == "m-dependent") {
    Mv <- attr(corstr, "Mv")
    
    #check that M is not too large
    if(Mv >= max(len)){
      stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
    }
    
    # Build block diagonal similar to in exchangeable case, also get row indices and column
    # indices for fast matrix updating later.
    tmp <- get_block_diag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    
  } else if(corstr == "unstructured"){
    
    if( max(len^2 - len)/2 > length(len)){
      stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
    }
    tmp <- get_block_diag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    
  } else if(corstr == "fixed") {
    corr.mat <- attr(corstr, "corr.mat")
    
    # check if matrix meets some basic conditions
    corr.mat <- check_fixed_mat(corr.mat, len)
    
    R.alpha.inv <- as(get_alpha_inv_fixed(corr.mat, len), "symmetricMatrix")/phi
    BlockDiag <- get_block_diag(len)$BDiag
    
  } else if(corstr == "userdefined") {
    corr.mat <- attr(corstr, "corr.mat")
  
    corr.mat <- check_user_mat(corr.mat, len)
    
    # get the structure of the correlation matrix in a way that
    # I can use later on.
    tmp1 <- getUserStructure(corr.mat)
    corr.list <- tmp1$corr.list
    user.row <- tmp1$row.vec
    user.col <- tmp1$col.vec
    struct.vec <- tmp1$struct.vec
    
    # the same block diagonal trick.
    tmp2 <- get_block_diag(len)
    BlockDiag <- tmp2$BDiag
    row.vec <- tmp2$row.vec
    col.vec <- tmp2$col.vec
    
  }
  
  
  
  ###############################################################################################
  # Fit model ###################################################################################
  ###############################################################################################
  
  done <- FALSE 
  converged <- FALSE
  count <- 1
  beta.old <- beta
  unstable <- FALSE
  phi.old <- phi
  
  # Fisher scoring loop
  while(!done && count <= maxit) {
 
    eta <- as.vector(x %*% beta) + offset
    mu <- InvLink(eta)
    diag(StdErr) <- sqrt(1/VarFun(mu))
    
    if(!scale.fix){
      phi <- update_phi(y, mu, VarFun, p, StdErr, included, includedlen, 
                       sqrtW, useP)
    }
    
    #Compute residuals
    Resid <- residuals_geelm(StdErr = StdErr, included = included, 
                            sqrtW = sqrtW, YY = y, mu = mu)
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(corstr == "ar1") { 
      
      alpha.new <- update_alpha_ar(y, mu, VarFun, phi, id, len, StdErr, Resid, p,
                                 included, includedlen, includedvec, allobs,
                                 sqrtW, BlockDiag, useP)
      R.alpha.inv <- get_alpha_inv_ar(alpha.new, a1, a2, a3, a4, row.vec, col.vec)/phi
      
    } else if(corstr == "exchangeable") {
     
      alpha.new <- update_alpha_ex(y, mu, VarFun, phi, id, len, StdErr,
                                 Resid, p, BlockDiag, included,
                                 includedlen, sqrtW, useP)
      R.alpha.inv <- get_alpha_inv_ex(alpha.new, n.vec, BlockDiag)/phi
      
    } else if(corstr == "m-dependent") {
      
      if(Mv==1){ #???IS THIS RIGHT???? - this should be done way earlier if Mv = 1 means AR1 (but is that right???)
        alpha.new <- update_alpha_ar(y, mu, VarFun, phi, id, len, StdErr, Resid, p,
                                   included, includedlen, includedvec, allobs,
                                   sqrtW, BlockDiag, useP)
      }else{
        alpha.new <- update_alpha_mdep(y, mu, VarFun, phi, id, len,
                                     StdErr, Resid, p, BlockDiag, Mv,
                                     included, includedlen,
                                     allobs, sqrtW, useP)
        if(sum(len>Mv) <= p){
          unstable <- TRUE
        }
      }
      if(any(alpha.new >= 1)){ 
        done <- TRUE
        warning("An estimated correlation great than 1 was found, stopping before convergence.")
      }
      R.alpha.inv <- get_alpha_inv_mdep(alpha.new, len, row.vec, col.vec)/phi
      
    } else if(corstr == "unstructured") {
     
      alpha.new <- update_alpha_unstruc(y, mu, VarFun, phi, id, len,
                                      StdErr, Resid,  p, BlockDiag,
                                      included, includedlen, allobs,
                                      sqrtW, useP)
      # This has happened to me (greater than 1 correlation estimate)
      if(any(alpha.new >= 1)){
        done <- TRUE
        warning("An estimated correlation great than 1 was found, stopping before convergence.")
      }
      R.alpha.inv <- get_alpha_inv_unstruc(alpha.new, len, row.vec, col.vec)/phi

    } else if(corstr == "fixed") {
      
      # FIXED CORRELATION, DON'T NEED TO RECOMPUTE
      R.alpha.inv <- R.alpha.inv*phi.old/phi 
      alpha.new <- NULL
     
    } else if(corstr == "userdefined") {
     
      alpha.new <- update_alpha_user(y, mu, phi, id, len, StdErr, Resid,
                                   p, BlockDiag, user.row, user.col,
                                   corr.list, included, includedlen,
                                   allobs, sqrtW, useP)
      R.alpha.inv <- get_alpha_inv_user(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
      
    } else if(corstr == "independence") {
      
      R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
      alpha.new <- "independent"
    }
    
    beta.list <- update_beta(y, x, beta, offset, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, 
                            dInvLinkdEta, tol, W, included)
    beta <- beta.list$beta
    
    phi.old <- phi 
    if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){
      converged <- TRUE
      done <- TRUE
    }
    beta.old <- beta
    
    count <- count + 1
  }
  
  biggest <- which.max(len)[1]
  index <- sum(len[1:biggest])-len[biggest]
  
  # Information for geem-class output
  if(K == 1){
    biggest.R.alpha.inv <- R.alpha.inv
    if(corstr == "fixed") {
      biggest.R.alpha <- corr.mat*phi
    } else {
      biggest.R.alpha <- solve(R.alpha.inv)
    }
  } else { # case: K != 1
    biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
    if(corstr == "fixed"){ 
      biggest.R.alpha <- corr.mat[1:len[biggest] , 1:len[biggest]]*phi
    }else{
      biggest.R.alpha <- solve(biggest.R.alpha.inv)
    }
  }
  
  eta <- as.vector(x %*% beta) + offset
  if(sandwich){
    sandvar.list <- get_sandwich(y, x, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, 
                                beta.list$hess, StdErr, dInvLinkdEta, BlockDiag, W, included)
  } else {
    sandvar.list <- list(sandvar = "no sandwich")
  }
  
  alpha <- alpha.new
  if(is.character(alpha)){
    alpha <- 0
  }
  if(corstr == "fixed") {
    alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat, 1) != 0)])
  }
  

  #Extra computations needed for geepack output 
  pearson_resid <- (y - mu) * diag(sqrtW)/VarFun(mu)
  fitted.values <- InvLink(eta)
  
  #Issue warnings
  if (!converged) {
    warning("Did not converge")
  }
  if (unstable) {
    warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")
  } 
  
  # Prep betas for output - add parameter names
  beta <- as.vector(beta)
  names(beta) <- colnames(x)
  
  # Collect and return everything
  out <- list(alpha = alpha,
              beta = beta,
              phi = phi,
              niter = count - 1,
              naiv.var = solve(beta.list$hess), 
              corr = corstr, 
              clusz = len,
              FunList = family,
              offset = offset,
              eta = eta,
              weights = weights,
              biggest.R.alpha = biggest.R.alpha/phi,
              resid = pearson_resid,
              fitted.values = fitted.values,
              vbeta =  as.matrix(sandvar.list$sandvar))
  class(out) <- c("geelm.fit", "list")
  
  out
  
  
}







#####################################################################################
## Not exported below
#####################################################################################


## Compute residuals
residuals_geelm <- function(StdErr, included, sqrtW, YY, mu) {
  StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
}