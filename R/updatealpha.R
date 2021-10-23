### Update the alpha (possibly) vector for the USERDEFINED correlation matrix.
update_alpha_user <- function(YY, mu, phi, id, len, StdErr, Resid,
                            p, BlockDiag, row.vec, col.vec, corr.list,
                            included, includedlen, allobs,sqrtW, useP){
  ml <- max(len)

  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", length(corr.list))

  index <- cumsum(len)-len

  for(i in 1:length(alpha.new)){
    newrow <- NULL
    newcol <- NULL
    for(j in 1:length(corr.list[[i]])){
      newrow <- c(newrow, index[which(len >= col.vec[corr.list[[i]]][j])] + row.vec[corr.list[[i]][j]])
      newcol <- c(newcol, index[which(len >= col.vec[corr.list[[i]]][j])] + col.vec[corr.list[[i]][j]])
    }

    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- phi*(length(newrow) - useP * p)
    }else{
      denom <- phi*(sum(bdtmp!=0) - useP * p)
    }
    alpha.new[i] <- sum(bdtmp)/denom
  }
  return(alpha.new)
}



### Calculate the parameter for the EXCHANGEABLE correlation structure
update_alpha_ex <- function(YY, mu, VarFun, phi, id, len, StdErr,
                          Resid, p, BlockDiag, included,
                          includedlen, sqrtW, useP){
  W <- sqrtW^2
  Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - mu)
  #print(mean(StdErr))
  #print(mean(YY-mu))
  
  
  #denom <- phi*(sum(triu(included %*% BlockDiag %*% included, k=1)) - useP * p)
  denom <- phi*(sum(includedlen*(includedlen-1))/2 - useP * p)
  
  #BlockDiag <- StdErr  %*%Diagonal(x = YY - mu) %*%  BlockDiag %*% W %*% included %*% Diagonal(x = YY - mu)  %*% StdErr
  BlockDiag <- StdErr  %*%Diagonal(x = YY - mu) %*%  included %*% BlockDiag %*% W %*% Diagonal(x = YY - mu)  %*% StdErr
    
  alpha <- sum(triu(BlockDiag, k=1))

  numpos <- sum(triu(BlockDiag,k=1) != 0)
  #print(numpos/2)
  #print(sum(triu(included %*% BlockDiag %*% included, k=1)))
  #print(sum(includedlen*(includedlen-1))/2)
  
  #alpha <- sum(triu(Resid %*% BlockDiag %*% Resid, k=1))
  alpha.new <- alpha/denom
  
  return(alpha.new)
}

### Calculate the parameters for the M-DEPENDENT correlation structure
update_alpha_mdep <- function(YY, mu, VarFun, phi, id, len,
                            StdErr, Resid, p, BlockDiag, m,
                            included, includedlen, allobs, sqrtW, useP){

  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", m)
  for(i in 1:m){
    if(sum(includedlen>i) > p){
      bandmat <- drop0(band(BlockDiag, i,i))

      if(allobs){
        alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i)) - useP * p))
      }else{
        alpha.new[i] <- sum( bandmat)/(phi*(length(bandmat@i) - useP * p))
      }

    }else{
      # If we don't have many observations for a certain parameter, don't divide by p
      # ensures we don't have NaN errors.
      bandmat <- drop0(band(BlockDiag, i,i))

      if(allobs){alpha.new[i] <- sum(bandmat)/(phi*(sum(as.numeric(len>i)*(len-i))))
      }else{
        alpha.new[i] <- sum( bandmat)/(phi*length(bandmat@i))
      }

    }
  }
  return(alpha.new)
}

### Calculate the parameter for the AR-1 correlation, also used for 1-DEPENDENT
update_alpha_ar <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid, p,
                          included, includedlen, includedvec, allobs,
                          sqrtW, BlockDiag, useP){

  Resid <- StdErr %*%  included %*% sqrtW %*% Diagonal(x = YY - mu)

  denom <- phi*(sum(band(triu(included %*% BlockDiag %*% included, k=1), k1=1, k2=1)) - useP * p)

  num <- sum(band(triu(Resid %*% BlockDiag %*% Resid, k=1), k1=1, k2=1))
  
  alpha <- num/denom
  return(as.numeric(alpha))
}


### Calculate alpha values for UNSTRUCTURED correlation
update_alpha_unstruc <- function(YY, mu, VarFun, phi, id, len, StdErr, Resid,
                               p, BlockDiag, included, includedlen,
                               allobs, sqrtW, useP){

  ml <- max(len)

  BlockDiag <- Resid %*% BlockDiag %*% Resid
  alpha.new <- vector("numeric", sum(1:(ml-1)))
  lalph <- length(alpha.new)

  row.vec <- NULL
  col.vec <- NULL
  for(i in 2:ml){
    row.vec <- c(row.vec, 1:(i-1))
    col.vec <- c(col.vec, rep(i, each=i-1))
  }
  index <- cumsum(len)-len
  if(sum(includedlen == max(len)) <= p){stop("Number of clusters of largest size is less than p.")}
  for(i in 1:lalph){
    # Get all of the indices of the matrix corresponding to the correlation
    # we want to estimate.
    newrow <- index[which(len>=col.vec[i])] + row.vec[i]
    newcol <- index[which(len>=col.vec[i])] + col.vec[i]
    bdtmp <- BlockDiag[cbind(newrow, newcol)]
    if(allobs){
      denom <- (phi*(length(newrow) - useP * p))
    }else{
      denom <- (phi*(sum(bdtmp!=0) - useP * p))
    }
    alpha.new[i] <- sum(bdtmp)/denom
  }

  return(alpha.new)
}



