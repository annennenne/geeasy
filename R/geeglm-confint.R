#Compute Wald confidence intervals for mean structure parameters of geeglm object

confint.geeglm <- function(object, parm = NULL, level = 0.95, std.err = "sandwich") {
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