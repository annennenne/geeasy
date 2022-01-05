#' @title Plot parameter estimates with 95 pct confidence intervals. 
#' 
#' @description Parameter estimates are plotted along with error bars
#'     indicating 95 pct confidence intervals.
#' 
#' @details
#' One or more models can be supplied, and if the parameters have the
#'     same names across models, they will be grouped together
#'     allowing for easy comparison.
#' 
#' Note that models can be given with or without names. If names are
#' supplied (see example below), these are printed in the plot
#' legend. Otherwise, the name of the model object is printed there
#' instead.
#' 
#' Implementation details: As a default, the estimates are extracted
#' from the \code{$coefficients} slot from the model object and
#' confidence intervals are computed by calling \code{confint()}. This
#' means that \code{plotEst} supports all models that have a
#' \code{coefficents} slot and a \code{confint} method.
#'   
#' @param ... One or more models. Currently, \code{lm}, \code{glm}, \code{geeglm},
#' \code{geelm}, and \code{mice} models are supported, but other model types may
#' be applicable as well, see details. 
#' 
#' @param intercept Logical indicating whether the intercept should be plotted (defaults 
#' to \code{TRUE}). 
#' 
#' @param par Which model parameters to plot estimates for, given as character strings. 
#' Default is \code{NULL} which means that all parameter estimates are plotted. 
#' 
#' @param colors Color scale to use if several models are plotted. Defaults to 
#' a color blind friendly scale. If there are more models than there are colors,
#' the color values are repeated. 
#'
#' @return No return values; called for side effects.
#' 
#' @author Anne Helby Petersen
#' 
#' @examples 
#' # Fit example models
#' data(iris)
#' m1 <- lm(Sepal.Length ~ Petal.Length, iris)
#' m2 <- lm(Sepal.Length ~ Petal.Length + Petal.Width, iris)
#' 
#' # Plot one model
#' plotEst(m2)
#' 
#' # Plot two models with default model labels
#' plotEst(m1, m2)
#' 
#' # Plot two models with custom model labels (simple)
#' plotEst(model1 = m1, model2 = m2)
#' 
#' # Plot two models with custom model labels (with spacing)
#' plotEst(`Simple model` = m1, `Full petal model` = m2)
#' 
#' # Plot two models without intercept
#' plotEst(m1, m2, intercept = FALSE)
#' 
#' # Plot two models with custom parameter subset
#' plotEst(m1, m2, par = c("Petal.Length"))
#' 
#' # Plot two models with custom color scale given by color names
#' plotEst(m1, m2, colors = c("red", "blue"))
#' 
#' # Plot two models with custom color scale given by color hex codes
#' #  note: only the first colors are used as there are more 
#' #  colors than models
#' plotEst(m1, m2, colors = c("#CC6666", "#9999CC", "#66CC99"))
#' 
#' 
#' @importFrom ggplot2 ggplot aes_string geom_point position_dodge coord_flip 
#'  scale_x_discrete ylab scale_color_manual geom_errorbar
#'
#' @importFrom stats confint
#' @export
plotEst <- function(..., intercept = TRUE, par = NULL, colors = NULL) {
  innames <- sapply(substitute(list(...))[-1], deparse)
  pF <- NULL
  ms <- list(...)
  mnames <- names(ms)
  
  if (is.null(mnames)) {
    mnames <- innames
  } else {
    emptyplaces <- mnames == ""
    mnames[emptyplaces] <- innames[emptyplaces]
  }
  
  n <- length(ms)
  
  for (i in 1:n) {
    m <- ms[[i]]  
    name <- mnames[[i]]
    if ("mipo" %in% class(m)) { 
      res <- summary(m) 
      if ("term" %in% names(res)) {
        sumci <- summary(m, conf.int = TRUE)
        mF <- data.frame(est = res$estimate,
                         lwr = sumci$`2.5 %`,
                         upr = sumci$`97.5 %`,
                         par = res$term,
                         model = name,
                         stringsAsFactors = FALSE)
      } else {
        #case: old version of mice, unclear when they made
        #this change
        mF <- data.frame(est = res[,1],
                         lwr = summary(m, conf.int = TRUE)[,6],
                         upr = summary(m, conf.int = TRUE)[,7],
                         par = rownames(res),
                         model = name,
                         stringsAsFactors = FALSE)
      }
    } else { #case: lm/glm/default
      mF <- data.frame(est = m$coefficients,
                       lwr = confint(m)[,1],
                       upr = confint(m)[,2],
                       par = names(coef(m)),
                       model = name,
                       stringsAsFactors = FALSE)
    }
    pF <- rbind(pF, mF)
  }
  
  pF$model <- factor(pF$model, levels = rev(mnames))
  
  if (!intercept) {
    pF <- pF[pF$par != "(Intercept)", ]
  }
  
  if (!is.null(par)) {
    pF <- pF[pF$par %in% par, ]
  }

  pF$par <- factor(pF$par)
  
  ## Handle colors
  ## Default: Color blind friendly palette
  ## Repeat colors if there are more models than there are colors
  if (is.null(colors)) {
    colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }
  nModels <- length(unique(pF$model))
  nCols <- length(colors)
  nrep <- ceiling(nModels/nCols) 
  colors <- rep(colors, nrep)
  

  q <- ggplot(pF, aes_string(x = "par", y = "est", 
                             ymin = "lwr", ymax = "upr", 
                             col = "model")) +
    geom_point(position = position_dodge(0.8)) + 
    geom_errorbar(position = position_dodge(0.8)) +
    coord_flip() + 
    scale_x_discrete("", limits = rev(levels(pF$par))) +
    ylab("Estimate") +
    scale_color_manual("Model", breaks = rev(levels(pF$model)),
                       values = colors)
                        
  q
}



