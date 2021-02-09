geem.fit <- function(x, y, id, offset, family, weights, control, corstr,
                     allobs, start = NULL) {
  
  # Unpack family functions
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
  
  # Number of covariates (p) and observations in total
  p <- dim(x)[2] 
  nn <- dim(x)[1] 
  
  ## Basic check to see if link and variance functions make any kind of sense
  meanY <- mean(y)
  linkOfMean <- LinkFun(meanY) - mean(offset)
  if(!is.finite(linkOfMean)) {
    stop("Infinite or NaN in the link of the mean of responses. Make sure link function makes sense for these data.")
  }
  if(!is.finite(VarFun(meanY))) {
    stop("Infinite or NaN in the variance of the mean of responses. Make sure variance function makes sense for these data.")
  }
  
  
  #!!!! quick solution, consider if this should be done more elegantly (changing code below)
  if (control$std.err == "san.se") {
    sandwich <- TRUE
  } else {
    sandwich <- FALSE
  }
  
  ###############################################################################################
  # Initialization
  ###############################################################################################
  
  #OLD BETA INITIALIZATION CODE - DELETE???? 
  ##############
  if (FALSE) {
    # Initialize beta
    beta <- control$init.beta
    
    if(is.null(beta)) {
      # Is there an intercept column?
      interceptcol <- apply(x == 1, 2, all)
      
      if(any(interceptcol)){
        #if there is an intercept and no initial beta, then use link of mean of response
        beta <- rep(0, p)
        beta[which(interceptcol)] <- linkOfMean
      } else {
        stop("Must supply an initial beta if not using an intercept.")
      }
    }
  }
  ###############
  
  #Initialize beta
  #First choice: Values supplied via control
  #Second choice: Values from glm fit
  beta <- control$init.beta
  if (is.null(beta)) {
    m_glm <- glm.fit(x = x, y = y, offset = offset, family = family, weights = weights)
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
  
  ###############################################################################################
  # Perform computations needed no matter which corstr is used
  ###############################################################################################
  
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
  Resid <- Diagonal(nn) #?? this is called from step below but never modified. 
                        #maybe move this definition into the estimator
  
  
  
  
  ###############################################################################################
  # Initialization according to corstr
  ###############################################################################################
  
  # Initialize for each correlation structure
  if(corstr == "independence") {
    R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
    BlockDiag <- getBlockDiag(len)$BDiag
  } else if(corstr == "ar1"){
    tmp <- buildAlphaInvAR(len)
    # These are the vectors needed to update the inverse correlation
    a1 <- tmp$a1
    a2 <- tmp$a2
    a3 <- tmp$a3
    a4 <- tmp$a4
    # row.vec and col.vec for the big block diagonal of correlation inverses
    # both are vectors of indices that facilitate in updating R.alpha.inv
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    BlockDiag <- getBlockDiag(len)$BDiag
    
  } else if(corstr == "exchangeable"){
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
    
  } else if(corstr == "m-dependent") {
    
    Mv <- attr(corstr, "Mv")
    
    #check that M is not too large
    if(Mv >= max(len)){
      stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
    }
    
    # Build block diagonal similar to in exchangeable case, also get row indices and column
    # indices for fast matrix updating later.
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    
  } else if(corstr == "unstructured"){
    
    if( max(len^2 - len)/2 > length(len)){
      stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
    }
    tmp <- getBlockDiag(len)
    BlockDiag <- tmp$BDiag
    row.vec <- tmp$row.vec
    col.vec <- tmp$col.vec
    
  } else if(corstr == "fixed") {
    corr.mat <- attr(corstr, "corr.mat")
    
    # check if matrix meets some basic conditions
    corr.mat <- checkFixedMat(corr.mat, len)
    
    R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
    BlockDiag <- getBlockDiag(len)$BDiag
    
  } else if(corstr == "userdefined") {
    corr.mat <- attr(corstr, "corr.mat")
  
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
    
  }
  
  
  #########################################################################
  # Fit model #############################################################
  #########################################################################

  
  
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
      phi <- updatePhi(y, mu, VarFun, p, StdErr, included, includedlen, 
                       sqrtW, useP)
    }
    
    #Compute residuals
    Resid <- residuals_geem(StdErr = StdErr, included = included, 
                            sqrtW = sqrtW, YY = y, mu = mu)
    
    
    ## Calculate alpha, R(alpha)^(-1) / phi
    if(corstr == "ar1") { 
      
      alpha.new <- updateAlphaAR(y, mu, VarFun, phi, id, len, StdErr, Resid, p,
                                 included, includedlen, includedvec, allobs,
                                 sqrtW, BlockDiag, useP)
      R.alpha.inv <- getAlphaInvAR(alpha.new, a1, a2, a3, a4, row.vec, col.vec)/phi
      
    } else if(corstr == "exchangeable") {
     
      alpha.new <- updateAlphaEX(y, mu, VarFun, phi, id, len, StdErr,
                                 Resid, p, BlockDiag, included,
                                 includedlen, sqrtW, useP)
      R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
      
    } else if(corstr == "m-dependent") {
      
      if(Mv==1){ #???? change - this should be done way earlier if Mv = 1 means AR1 (but is that right???)
        alpha.new <- updateAlphaAR(y, mu, VarFun, phi, id, len, StdErr, Resid, p,
                                   included, includedlen, includedvec, allobs,
                                   sqrtW, BlockDiag, useP)
      }else{
        alpha.new <- updateAlphaMDEP(y, mu, VarFun, phi, id, len,
                                     StdErr, Resid, p, BlockDiag, Mv,
                                     included, includedlen,
                                     allobs, sqrtW, useP)
        if(sum(len>Mv) <= p){
          unstable <- TRUE
        }
      }
      if(any(alpha.new >= 1)){ #!!! this check should probably be done in all cases?
        done <- TRUE
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, col.vec)/phi
      
    } else if(corstr == "unstructured") {
     
      alpha.new <- updateAlphaUnstruc(y, mu, VarFun, phi, id, len,
                                      StdErr, Resid,  p, BlockDiag,
                                      included, includedlen, allobs,
                                      sqrtW, useP)
      # This has happened to me (greater than 1 correlation estimate)
      if(any(alpha.new >= 1)){
        done <- TRUE
        warning("some estimated correlation is greater than 1, stopping.")
      }
      R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, row.vec, col.vec)/phi

    } else if(corstr == "fixed") {
      
      # FIXED CORRELATION, DON'T NEED TO RECOMPUTE
      #!! Note: if !scale.fix phi may be updated. But otherwise this is not happening. 
      #!! Check if we can avoid iterative estimation here altogether and whether phi
      #   phi is actually updated when corstr is "fixed". 
      R.alpha.inv <- R.alpha.inv*phi.old/phi 
      alpha.new <- NULL
     
    } else if(corstr == "userdefined") {
     
      alpha.new <- updateAlphaUser(y, mu, phi, id, len, StdErr, Resid,
                                   p, BlockDiag, user.row, user.col,
                                   corr.list, included, includedlen,
                                   allobs, sqrtW, useP)
      R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
      
    } else if(corstr == "independence") {
      ###!!! use glm in this scenario? don't go into recursive updating here? 
      R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
      alpha.new <- "independent"
    }
    
    
    beta.list <- updateBeta(y, x, beta, offset, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, 
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
    sandvar.list <- getSandwich(y, x, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, 
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
  } #!!! consider restructuring so that this slot is used more generally for error-messages,
  # which would allow for it to be reused by other correlation structures
  
  # Prep betas for output - add parameter names
  beta <- as.vector(beta)
  names(beta) <- colnames(x)
  
  # Collect and return everything
  out <- list(alpha = alpha,
              beta = beta,
              phi = phi,
              niter = count - 1,
             # converged = converged,
            #  unstable = unstable,
              naiv.var = solve(beta.list$hess), 
             # var = sandvar.list$sandvar,
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
  
  out
  
  
}




## Compute residuals
residuals_geem <- function(StdErr, included, sqrtW, YY, mu) {
  StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
}