identicalLists <- function(l1, l2, ignore = NULL) {
  len1 <- length(l1)
  len2 <- length(l2)
  
  if(len1 != len2) {
    message("Different list lengths")
    return(FALSE)
  }
  
  if (!is.null(ignore)) {
    l1[ignore] <- NULL
    l2[ignore] <- NULL
    len1 <- length(l1)
  }
  out <- rep(NA, len1)
  for (i in 1:len1) {
    out[i] <- identical(l1[[i]], l2[[i]])
  }
  if (all(out)) {
    return(TRUE)
  } else {
    message(paste("problems at:", paste(names(l1)[which(!out)], collapse = ", ")))
    return(FALSE)
  }
  out
}
