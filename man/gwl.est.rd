\name{gwl.est}
\alias{gwl.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographically weighted lasso
}
\description{
This function fits a geographically weighted lasso (GWL) model
}
\usage{
gwl.est(form, locs, data, kernel = "exp", cv.tol)
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
This function estimates penalized spatially varying coefficients using the geographically weighed regression
and lasso approaches. Spatial kernel weights are applied to observations using the estimated kernel bandwidth to 
estimate local models at each data point. The kernel bandwidth and lasso solutions are currently estimated using 
cross-validation with an exponential or Gaussian kernel function. Some regression coefficients may be penalized to zero. 
The function estimates regression coefficients, the outcome variable values, and the model fit.
}
\value{
 A list with the following items:

  \item{phi }{Kernel bandwidth}
  \item{RMSPE }{Root mean squared prediction error from bandwidth estimation}
  \item{beta }{Matrix of estimated regression coefficients, where a row contains the coefficients
  for one regression term for all data points}
  \item{yhat }{Estimated outcome variable values}
  \item{RMSE }{Root mean squared error from estimation}
  \item{rsquare }{Approximate R-square for GWR model}
  
}
\references{
Wheeler DC (2009) Simultaneous coefficient penalization and model selection in geographically weighted regression: 
The geographically weighted lasso. Environment and Planning A, 41: 722-742
}
\author{
David Wheeler
}


\seealso{
  \code{\link{gwrr.est}}
}
\examples{
data(columbus)
locs <- cbind(columbus$x, columbus$y)
col.gwl <- gwl.est(crime ~ income + houseval, locs, columbus, "exp")
plot(col.gwl$beta[2,], col.gwl$beta[3,])
plot(columbus$x, columbus$y, cex=col.gwl$beta[1,]/10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

