#' Control estimation of GEE models
#' 
#' Settings for controlling technical details of GEE fitting via geelm. 
#' 
#' @inheritParams geepack::geese.control
#' @inheritParams geeM::geem
#' 
#' @param tol Tolerance for asserting convergence. 
#' 
#' @param useP If set to \code{FALSE}, do not use the n-p correction for 
#'    dispersion and correlation estimates. This can be 
#'    useful when the number of observations is small, as subtracting p may yield 
#'    correlations greater than 1.
#'    
#' @param std.err Character string specifying which standard error
#'     estimation method should be used. Supported options are
#'     `san.se` (sandwich SE) and `naive`.
#'
#' @return A list of values used for controlling model fitting.
#' 
#' @export
geelm.control <- function(#init.alpha = NULL, 
                         init.beta = NULL, init.phi = 1,
                         tol = 0.00001, maxit = 20, scale.fix = FALSE,
                         useP = TRUE, std.err = "san.se") {
  
  
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix = TRUE, then init.phi must be supplied")
  }
  
  useP <- as.numeric(useP)
  
  list(init.beta = init.beta, 
       init.phi = init.phi,
       tol = tol,
       maxit = maxit,
       scale.fix = scale.fix,
       useP = useP,
       jack = 0, #for geese methods
       j1s = 0, #for geese methods
       fij = 0, #for geese methods
       std.err = std.err
       )
}






