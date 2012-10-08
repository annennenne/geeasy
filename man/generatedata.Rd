\name{generatedata}
\alias{generatedata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates simulated count data
}
\description{
This function returns a data frame that can be used by \code{fitGEE}
}
\usage{
generatedata(beta, alpha, gamma, X, T, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
a numeric vector
}
  \item{alpha}{
a numeric scalar
}
  \item{gamma}{
a numeric scalar
}
  \item{X}{
a numeric design matrix
}
  \item{T}{
a numeric scalar, the number of observations per subject.
}
  \item{n}{
a numeric scalar, the number of subjects.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a data frame.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Nicholas Henderson
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
X <- cbind(rep(1,5),c(-.5,-.25,0,.25,.5))
testdat <- generatedata(beta=c(1,.5),alpha=.2,gamma=.5,X=X,T=5,n=3000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
