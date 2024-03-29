% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drop1.R
\name{drop1.geem}
\alias{drop1.geem}
\title{Drop All Possible Single Terms to a \code{geem} Model Using Wald or Score Test}
\usage{
\method{drop1}{geem}(
  object,
  scope,
  test = c("Wald", "none", "score", "sasscore"),
  method = c("robust", "naive", "sandwich"),
  ...
)
}
\arguments{
\item{object}{a fitted object of class geese.}

\item{scope}{a formula giving the terms to be considered for adding
or dropping.}

\item{test}{the type of test to include.}

\item{method}{Indicates which method is used for computing the
standard error. \code{robust} is the default and corresponds to
the modified sandwich estimator. \code{naive} is the classical
naive variance estimate.  \code{sandwich} is an alias for
\code{robust}.}

\item{\dots}{other arguments. Not currently used}
}
\value{
An object of class \code{anova} summarizing the differences in
fit between the models.
}
\description{
Compute all the single terms in the scope argument that can dropped from the
model, and compute a table of the corresponding Wald test statistics.
}
\examples{

library(geeM)
library(geepack)
data(ohio)
fit <- geem(resp ~ age + smoke + age:smoke, id=id, data=ohio,
            family="binomial", corstr="exch", scale.fix=TRUE)
drop1(fit)

}
\seealso{
\code{\link{drop1}}, \code{geem}
}
\author{
Claus Ekstrom \email{claus@ekstroem.dk}
}
\keyword{htest}
