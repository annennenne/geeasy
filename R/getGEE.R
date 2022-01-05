#' Get information for a geelm/geeglm object
#' 
#' @param object A geeglm model object as obtained from \code{geepack::geeglm()}
#' or a geelm object as obtained from \code{geeasy::geelm()}
#' 
#' @param name Name of the slot/component of the geelm/geeglm object that should be returned.
#' See list of possible names in details below.
#' 
#' @param ... Any additional arguments passed on to other functions.
#' 
#' @details The allowed names are:
#' 
#' \code{coefficients} or \code{beta}: Coefficients from the mean
#' structure model (betas) on their original scales
#' 
#' \code{residuals}: Pearson residuals, in the order of the inputted
#' dataset (with NAs omitted).
#' 
#' \code{fitted.values}: Fitted values (response scale), in the order
#'  of the inputted dataset (with NAs omitted).
#'  
#' \code{rank}: The rank of the model matrix, i.e. the number of
#' estimated mean structure coefficients.
#' 
#' \code{family}: A family object specifying which exponential family
#' was used for fitting the mean structure model, see
#' \code{\link{family}} for more information.
#' 
#' \code{linear.predictors}: The linear predictor on the original
#' scale.
#' 
#' \code{df.residuals}: Residual degrees of freedom.
#' 
#' \code{corstr}: Name of the correlation structured imposed on the
#' model. If the correlation structure requires further information,
#' it is stored in a suitably named attribute. For example, for
#' m-dependent correlation structures, the m scalar is available in an
#' attribute named \code{Mv}.
#' 
#' \code{std.err}: Method used to estimate the standard error of the
#' mean structure coefficients (betas).
#' 
#' \code{alpha}: The estimated parameter(s) for the variance
#' structure.
#' 
#' \code{gamma} or \code{dispersion}: The estimated dispersion
#' parameter.
#' 
#' \code{vbeta}: The estimated variance matrix of the mean structure
#' (beta) coefficients.
#' 
#' \code{beta.se}: The standard errors of the mean structure (beta)
#' coefficients.
#' 
#' \code{clusz} or \code{clustersizes}: Sizes of each of the clusters
#' in the data.
#' 
#' \code{nclusters}: The total number of clusters.
#'
#' @return The requested slot of the object.
#' 
#' @export
getGEE <- function(object, name, ...) {
  
  # All supported names that are slots in geelm/geeglm object
    geelm_slots <- c("coefficients", "beta",
                     "residuals", 
                     "fitted.values", 
                     "rank", 
                     "family", 
                     "linear.predictors",
                     "df.residuals",
                     "corstr", 
                     "std.err")
    
    ## All supported names that are slots in geese object 
    geese_slots <- c("alpha", 
                     "gamma", "dispersion",
                     "vbeta", 
                     "clusz", "clustersizes")
    
    ## All supported names that are not slots, but must be computed
    other_names <- c("beta.se", "nclusters")
  
    ## Store main object class and throw error if it's not geelm or geeglm
    objectclass <- class(object)
    if ("geelm" %in% objectclass) {
        objectclass <- "geelm"
    } else if ("geeglm" %in% objectclass) {
        objectclass <- "geeglm"
    } else stop("Input object must be of class either geelm or geeglm.")
    
    ## Check if name is supported
    if (!(name %in% c(geelm_slots, geese_slots, other_names))) {
        stop(paste(name, "is not the name of a geeglm component."))
    }
  
    ## Find geelm slot
    if (name %in% geelm_slots) {
        if (name == "beta") name <- "coefficients"
        return(object[[name]])
    } 
  
    ## Find geese slot
    if (name %in% geese_slots) {
        if (name == "dispersion") name <- "gamma"
        if (name == "clustersizes") name <- "clusz"    
        return(object$geese[[name]])
    }
    
    ## Compute other names
    if (name == "beta.se") {
        return(sqrt(diag(object$geese$vbeta)))
    }
    
    if (name == "nclusters") {
        return(length(object$geese$clusz))
    }
    
}
