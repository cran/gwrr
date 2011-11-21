\name{gwr.vdp}
\alias{gwr.vdp}
\title{
Collinearity diagnostics for geographically weighted regression
}
\description{
Uses the collinearity diagnostic tools of variance-decomposition proportions and condition indexes for
geographically weighted regression (GWR) models.
}
\usage{
gwr.vdp(form, locs, data, phi, kernel = "exp", sel.ci = 30, sel.vdp = 0.5)
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
  \item{phi}{
  The kernel bandwidth used in the GWR model
}
  \item{kernel}{
  The kernel weighting function used in the GWR model, either exp or gauss; exp is the default
}
  \item{sel.ci}{
  The threshold value to use for the condition index to indicate observations with a collinearity issue;
  indexes above this value will be flagged; the default is 30
}
  \item{sel.vdp}{
  The threshold value to use for the variance-decomposition proportion to indicate observations with a collinearity issue;
  proportions above this value will be flagged; the default is 0.5
}
}
\details{
This function calculates the variance-decomposition proportions and the condition indexes for the 
weighted design matrix used in a GWR model. The kernel function and bandwidth used to estimate the 
GWR model must be input to this function. Observations with a large condition index and relatively
large variance-decomposition proportions for more than one regression term indicate an issue with 
collinearity.
}
\value{
 A list with the following items:

  \item{condition }{Largest condition index for each observation}
  \item{vdp }{Variance-decomposition proportions for the largest variance component for each observation}
  \item{flag.cond }{True if largest condition index exceeds threshold}
  \item{flag.vdp }{True if variance-decomposition proportions for more than one term exceed threshold}
  \item{flag.cond.vdp }{True if condition index and variance-decompostion proportions exceed thresholds}
   
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
col.bw <- gwr.bw.est(crime ~ income + houseval, locs, columbus, "exp")
col.vdp <- gwr.vdp(crime ~ income + houseval, locs, columbus, col.bw$phi, "exp")
hist(col.vdp$condition)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gwr }
