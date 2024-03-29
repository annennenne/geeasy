% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotEst.R
\name{plotEst}
\alias{plotEst}
\title{Plot parameter estimates with 95 pct confidence intervals.}
\usage{
plotEst(..., intercept = TRUE, par = NULL, colors = NULL)
}
\arguments{
\item{...}{One or more models. Currently, \code{lm}, \code{glm}, \code{geeglm},
\code{geelm}, and \code{mice} models are supported, but other model types may
be applicable as well, see details.}

\item{intercept}{Logical indicating whether the intercept should be plotted (defaults
to \code{TRUE}).}

\item{par}{Which model parameters to plot estimates for, given as character strings.
Default is \code{NULL} which means that all parameter estimates are plotted.}

\item{colors}{Color scale to use if several models are plotted. Defaults to
a color blind friendly scale. If there are more models than there are colors,
the color values are repeated.}
}
\value{
No return values; called for side effects.
}
\description{
Parameter estimates are plotted along with error bars
indicating 95 pct confidence intervals.
}
\details{
One or more models can be supplied, and if the parameters have the
same names across models, they will be grouped together
allowing for easy comparison.

Note that models can be given with or without names. If names are
supplied (see example below), these are printed in the plot
legend. Otherwise, the name of the model object is printed there
instead.

Implementation details: As a default, the estimates are extracted
from the \code{$coefficients} slot from the model object and
confidence intervals are computed by calling \code{confint()}. This
means that \code{plotEst} supports all models that have a
\code{coefficents} slot and a \code{confint} method.
}
\examples{
# Fit example models
data(iris)
m1 <- lm(Sepal.Length ~ Petal.Length, iris)
m2 <- lm(Sepal.Length ~ Petal.Length + Petal.Width, iris)

# Plot one model
plotEst(m2)

# Plot two models with default model labels
plotEst(m1, m2)

# Plot two models with custom model labels (simple)
plotEst(model1 = m1, model2 = m2)

# Plot two models with custom model labels (with spacing)
plotEst(`Simple model` = m1, `Full petal model` = m2)

# Plot two models without intercept
plotEst(m1, m2, intercept = FALSE)

# Plot two models with custom parameter subset
plotEst(m1, m2, par = c("Petal.Length"))

# Plot two models with custom color scale given by color names
plotEst(m1, m2, colors = c("red", "blue"))

# Plot two models with custom color scale given by color hex codes
#  note: only the first colors are used as there are more 
#  colors than models
plotEst(m1, m2, colors = c("#CC6666", "#9999CC", "#66CC99"))


}
\author{
Anne Helby Petersen
}
