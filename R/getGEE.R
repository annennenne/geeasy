
#' Get information for a geelm/geeglm object
#' 
#' @param object A geeglm model object as obtained from \code{geepack::geeglm()}
#' or a geelm object as obtained from \code{geek::geelm()}
#' 
#' @param name Name of the slot/component of the geelm/geeglm object that should be returned.
#' See list of possible names in details below.
#' 
#' @details The allowed names are
#' 
#' * \code{coefficients}: Coefficients from the mean structure model (betas) on their 
#' original scales
#' * \code{residuals}:
#' * \code{fitted.values}: 
#' * \code{rank}: The rank of the model matrix, i.e. the number of estimated mean structure
#' coefficients.
#' * \code{qr}:
#' * \code{family}: A family object specifying which exponential family was used for fitting
#' the mean structure model, see \code{\link{stats::family}} for more information. 
#' * \code{linear.predictors}:
#' * \code{weights}: 
#' * \code{prior.weights}: The original weights used to produce this geeglm object (set
#' by user or defaulted to 1 for all observations).
#' * \code{df.residuals}:
#' * \code{y}:
#' * \code{model}:
#' * \code{call}: The original function call that produced this geeglm object.
#' * \code{formula}: 
#' * \code{terms}:
#' * \code{data}: The original dataset that was used for producing this geeglm object.
#' * \code{offset}:
#' * \code{control}: 
#' * \code{method}:
#' * \code{contrasts}:
#' * \code{xlevels}:
#' * \code{geese}: An object containing further information about the variance estimation, 
#' including a variance matrix for the beta-coefficients (\code{$vbeta}), the estimated 
#' coefficients for the working correlation matrix (\code{$alpha}), the estimated dispersion 
#' parameter (\code{$gamma}), and the individual cluster sizes (\code{$clusz}). See 
#' \code{\link{geepack::geese}} for more information. 
#' * \code{modelInfo}:
#' * \code{id}: 
#' * \code{corstr}: Name of the correlation structured imposed on the model. If the 
#' correlation structure requires further information, it is stored in a suitably named
#' attribute. For example, for m-dependent correlation structures, the m scalar is available
#' in an attribute named \code{Mv}. 
#' * \code{cor.link}:
#' * \code{std.err}: 
#' 
#' @export
getGEE <- function(object, name) {
  ALLSLOTS <- c("coefficients", "residuals", 
                "fitted.values", "rank", "qr", 
                "family", "linear.predictors", 
                "weights", "prior.weights", 
                "df.residuals", "y", "model", 
                "call", "formula", "terms", 
                "data", "offset", "control", 
                "method", "contrasts", "xlevels", 
                "geese", "modelInfo", "id", "corstr", 
                "cor.link", "std.err")
  if (!(name %in% ALLSLOTS)) {
    stop(paste(name, "is not a slot in a geeglm object."))
  }
  
  geeglm[[name]]
}