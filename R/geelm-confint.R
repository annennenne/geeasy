#' Confidence Intervals for geelm objects
#' 
#' Compute Wald confidence intervals for mean structure parameters of geelm object.
#' 
#' @param object a fitted model object.
#'
#' @param parm  specification of which parameters are to be given
#'           confidence intervals, either a vector of numbers or a vector
#'           of names.  If missing, all parameters are considered.
#'
#' @param level the confidence level required.
#'
#' @param std.err Which standard error estimation method that should
#'     be used for computing the confidence intervals. Only `san.se`
#'     is supported for geelm objects but `jack`, `j1s` or `fij` may
#'     be used for geeglm objects (if they have been estimated when
#'     fitting the model).
#' 
#' @param ... additional argument(s) for methods.
#' 
#' @return A matrix (or vector) with columns giving lower and upper
#'     confidence limits for each parameter.
#' 
#' @export
confint.geelm <- function(object, parm = NULL, level = 0.95, std.err = "san.se", ...) {
  betas <- object$coefficients
  if (std.err %in% c("san.se", "sandwich")) {
    v_betas <- object$geese$vbeta
  } else if (std.err == "jack") {
    v_betas <- object$geese$vbeta.ajs
  } else if (std.err == "j1s") {
    v_betas <- object$geese$vbeta.j1s
  } else if (std.err == "fij") {
    v_betas <- object$geese$vbeta.fij
  }
  se_betas <- sqrt(diag(v_betas))
  
  lwr_p <- (1 - level)/2
  upr_p <- 1 - lwr_p
  qn <- qnorm(upr_p)
  
  lwrs <- betas - qn * se_betas
  uprs <- betas + qn * se_betas
  
  out <- matrix(c(lwrs, uprs), nrow = length(betas), ncol = 2,
                dimnames = list(names(betas), 
                                paste(round(100 * c(lwr_p, upr_p), 1), " %", sep = "")))
  
  if (!is.null(parm)) {
    out <- out[parm, , drop = FALSE]
  }
  
  out  
}


#' @describeIn confint.geelm
#' 
#' @export
confint.geeglm <- function(object, parm = NULL, level = 0.95, std.err = "san.se", ...) {
  confint.geelm(object = object, parm = parm, level = level, std.err = std.err, ...)
}
