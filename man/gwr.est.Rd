\name{gwr.est}
\alias{gwr.est}
\title{
Geographically weighted regression
}
\description{
This function fits a geographically weighted regression (GWR) model
}
\usage{
gwr.est(form, coords, data, kernel = "exp", bw = TRUE, cv.tol)
}

\arguments{
  \item{form}{
  A regression model forumula, as in the functions lm and glm
}
  \item{coords}{
  A matrix of spatial coordinates of data points, where the x coordinate is first, then the y coordinate;
  coordinates are assumed to not be latitude and longitude, as Euclidean distance is calculated from coordinates
}
  \item{data}{
  A data frame with data to fit model
}
  \item{kernel}{
  A kernel weighting function, either exp or gauss, where exponential function is default
}
  \item{bw}{
  Either TRUE to estimate a bandwidth for the kernel function, or the bandwidth to use to fit the model;
  bandwidth is estimated by default
}
  \item{cv.tol}{
  A stopping tolerance in terms of cross-validation error for the bi-section search routine to estimate the kernel bandwidth using cross-validation;
  if missing an internally calculated value is used   
}
}
\details{
This function estimates spatially varying coefficients using the GWR approach.
Spatial kernel weights are applied to observations using the esimated or supplied kernel bandwidth to 
estimate local models at each data point. The bandwidth is currently estimated with cross-validation with
an exponential or Gaussian kernel function. The function estimates regression coefficients, the 
outcome variable values, and the model fit.
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
Wheeler DC (2007) Diagnostic tools and a remedial method for collinearity in geographically weighted regression.
Environment and Planning A, 39: 2464-2481.

Fotheringham AS, Brunsdon C, Charlton M (2002) Geographically weighted regression: The analysis of spatially varying relationships. 
John Wiley & Sons, West Sussex. 
}
\author{
David Wheeler
}

\seealso{
  \code{\link{gwr.bw.est}}
}
\examples{
data(columbus)
locs <- cbind(columbus$x, columbus$y)
col.gwr <- gwr.est(crime ~ income + houseval, locs, columbus, "exp")
plot(col.gwr$beta[2,], col.gwr$beta[3,])
plot(columbus$x, columbus$y, cex=col.gwr$beta[1,]/10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gwr }
