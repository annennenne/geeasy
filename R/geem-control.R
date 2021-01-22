geem.control <- function(#init.alpha = NULL, 
                         init.beta = NULL, init.phi = 1,
                         tol = 0.00001, maxit = 20, scale.fix = FALSE,
                         useP = TRUE) {
  
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix = TRUE, then init.phi must be supplied")
  }
  
  useP <- as.numeric(useP)
  
  list(init.beta = init.beta, 
       init.phi = init.phi,
       tol = tol,
       maxit = maxit,
       scale.fix = scale.fix,
       useP = useP)
}
























# Set the initial alpha value ##???? NONE OF THIS IS USED BELOW????
#*#  if(is.null(init.alpha)) {
#*#    alpha.new <- 0.2
#*#    #*# if(cor.match==4){
#*#    if(corstr == "m-dependent") {
#*#      # If corstr = "m-dep"
#*#      alpha.new <- 0.2^(1:Mv)
#*#      #*# }else if(cor.match==5) {
#*#    }else if(corstr == "unstructured") {
#*#      # If corstr = "unstructured"
#*#    alpha.new <- rep(0.2, sum(1:(max(len)-1)))
#*#    #*# } else if(cor.match==7) {
#*#  } else if(corstr == "userdefined") {
#*#    # If corstr = "userdefined"
#*#    alpha.new <- rep(0.2, max(unique(as.vector(corr.mat))))
#*#  }
#*#} else {
#*#  alpha.new <- init.alpha
#*#}

