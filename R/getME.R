#' @describeIn getGEE
#' 
#' @export
getME <- function(object, name, ...) {
  UseMethod("getME", object)
}

# Default: Use the original lme4 getME function
#'
#' @importFrom lme4 getME
#' @export
getME.default <- function(object, name, ...) {
  lme4::getME(object = object, name = name, ...)
}


#' 
#' @export
getME.geeglm <- function(object, name, ...) {
 getGEE(object = object, name = name) 
}


#' 
#' @export
getME.geelm <- function(object, name, ...) {
  getGEE(object = object, name = name)
}

