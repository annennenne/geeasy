getfam <- function(family){
  if(is.character(family)){
    family <- get(family, mode = "function", envir = parent.frame(2))  
  }
  
  if(is.function(family)){
    family <- family()
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family) && !is.null(family$family)){
    LinkFun <- family$linkfun
    InvLink <- family$linkinv
    VarFun <- family$variance
    InvLinkDeriv <- family$mu.eta
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv)    
  }else if(is.list(family)){
    if(length(match(names(family), c("LinkFun", "VarFun", "InvLink", "InvLinkDeriv"))) == 4){
      LinkFun <- family$LinkFun
      InvLink <- family$InvLink
      VarFun <- family$VarFun
      InvLinkDeriv <- family$InvLinkDeriv
    }else{
      LinkFun <- family[[1]]
      VarFun <- family[[2]]
      InvLink <- family[[3]]
      InvLinkDeriv <- family[[4]]
    }
    FunList <- list("LinkFun" = LinkFun, "VarFun" = VarFun, "InvLink" = InvLink, "InvLinkDeriv" = InvLinkDeriv) 
    return(FunList)
  }else{
    stop("problem with family argument: should be string, family object, or list of functions")
  }
}

### Get a block diagonal matrix. Each block has dimension corresponding to
### each cluster size.  By default, each block is just a matrix filled with ones.
getBlockDiag <- function(len, xvec=NULL){
  K <- length(len)
  if(is.null(xvec)){
    xvec <- rep.int(1, sum(len^2))
  }
  
  row.vec <- col.vec <- vector("numeric", sum(len^2))
  add.vec <- cumsum(len) - len
  index <- c(0, (cumsum(len^2) -len^2)[2:K], sum(len^2)) 
  for(i in 1:K){
    row.vec[(index[i] + 1):(index[i+1])] <- rep.int( (1:len[i]) + add.vec[i], len[i])
    col.vec[(index[i] + 1):(index[i+1])] <- rep( (1:len[i]) + add.vec[i], each=len[i])
  }	
  BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
  return(list(BDiag = as(BlockDiag, "symmetricMatrix"), row.vec =row.vec, col.vec=col.vec))
}


### Check some conditions on the FIXED correlation structure.
checkFixedMat <- function(corr.mat, len){
  if(is.null(corr.mat)){
    stop("corr.mat must be specified if using fixed correlation structure")
  }
  if(dim(corr.mat)[1] < max(len)){
    stop("Dimensions of corr.mat must be at least as large as largest cluster")
  }
  if(!isSymmetric(corr.mat)){
    stop("corr.mat must be symmetric")
  }
  if(determinant(corr.mat, logarithm=T)$modulus == -Inf){
    stop("supplied correlation matrix is not invertible.")
  }	
  return(corr.mat[1:max(len), 1:max(len)])	
}





### fitted function for geem object
fitted.geem <- function(object, ...){
  InvLink <- object$FunList$InvLink
  return(InvLink(object$eta))
}

predict.geem <- function(object, newdata = NULL,...){
  coefs <- object$beta
  if(is.null(newdata)){
    return(as.vector(object$X %*% object$beta))
  }else{
    if(dim(newdata)[2] != length(coefs)){warning("New observations must have the same number of rows as coefficients in the model")}
    return(as.vector(newdata %*% object$beta))
  }
}

coef.geem <- function(object, ...){
  coefs <- object$beta
  names(coefs) <- object$coefnames
  return(coefs)
}
