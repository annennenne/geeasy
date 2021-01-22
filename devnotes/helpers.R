identicalLists <- function(l1, l2, ignore = NULL, numIfNotExact = TRUE,
                           chatty = FALSE) {
  if (!is.null(ignore)) {
    l1[ignore] <- NULL
    l2[ignore] <- NULL
  }
  
  len1 <- length(l1)
  len2 <- length(l2)
  
  if(len1 != len2) {
    message("Different list lengths")
    return(FALSE)
  }
  
  out <- rep(NA, len1)
  for (i in 1:len1) {
    res <- identical(l1[[i]], l2[[i]])
    
    if (numIfNotExact & !res & is.numeric(l1[[i]])) {
      res <- max(abs(l1[[i]] - l2[[i]]))
      if (chatty) {
        print(paste("Diff for ", names(l1)[i], ": ", res, sep = ""))
      }
    } 
      
    out[i] <- as.numeric(res)
  }
  
  if (all(out == 1)) {
    return(TRUE)
  } else if (all(out == 0)) {
    message(paste("problems at:", paste(names(l1)[which(!out)], collapse = ", ")))
    return(FALSE)
  } else {
    return(max(out[out != 1]))
  }
}
