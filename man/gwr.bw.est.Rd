\name{gwr.bw.est}
\alias{gwr.bw.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation estimation of kernel bandwidth
}
\description{
Estimate the kernel function bandwidth with cross-validation
}
\usage{
gwr.bw.est(form, locs, data, kernel = "exp", cv.tol)
}

\arguments{
  \item{form}{
  A regression model forumula, as in the functions lm and glm
}
  \item{locs}{
  A matrix of spatial coordinates of data points, where the x coordinate is first, then the y coordinate;
  coordinates are assumed to not be latitude and longitude, as Euclidean distance is calculated from coordinates
}
  \item{data}{
  A data frame with data to fit model
}
  \item{kernel}{
  A kernel weighting function, either exp or gauss, where exponential function is default
}
  \item{cv.tol}{
  A stopping tolerance in terms of cross-validation error for the bi-section search routine to estimate the kernel bandwidth using cross-validation;
  if missing an internally calculated value is used   
}
}
\details{
This function estimates the kernel bandwidth in a GWR model with leave-one-out cross-validation.
It does not estimate the final regression coefficients or outcome variable.
}
\value{
 A list with the following items:

  \item{phi }{Kernel bandwidth}
  \item{RMSPE }{Root mean squared prediction error from bandwidth estimation}
  \item{cv.score }{Sum of squared prediction errors from bandwidth estimation}
}
\references{
Wheeler DC (2007) Diagnostic tools and a remedial method for collinearity in geographically weighted regression.
Environment and Planning A, 39: 2464-2481 
}
\author{
David Wheeler
}


\seealso{
  \code{\link{gwr.est}}
}
\examples{
data(columbus)
locs <- cbind(columbus$x, columbus$y)
col.bw <- gwr.bw.est(crime ~ income + houseval, locs, columbus, "exp")
col.gwr <- gwr.est(crime ~ income + houseval, locs, columbus, "exp", bw=col.bw$phi)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

