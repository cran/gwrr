\name{gwrr-package}
\alias{gwrr-package}
\alias{gwrr}
\docType{package}
\title{
Geographically weighted regression models with penalties and diagnostic tools
}
\description{
Fits geographically weighted regression (GWR) models and has tools to diagnose collinearity in the GWR models.
}
\details{
\tabular{ll}{
Package: \tab gwrr\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2011-11-21\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
David Wheeler

Maintainer: David Wheeler <dcwheels@gmail.com>

}
\references{
Wheeler DC (2007) Diagnostic tools and a remedial method for collinearity in geographically weighted regression.
Environment and Planning A, 39: 2464-2481. 

Fotheringham AS, Brunsdon C, Charlton M (2002) Geographically weighted regression: The analysis of spatially varying relationships. 
John Wiley & Sons, West Sussex.
}

\keyword{ package }

\examples{
data(columbus)
locs <- cbind(columbus$x, columbus$y)
col.gwr <- gwr.est(crime ~ income + houseval, locs, columbus, "exp")
plot(col.gwr$beta[2,], col.gwr$beta[3,])
col.vdp <- gwr.vdp(crime ~ income + houseval, locs, columbus, col.gwr$phi, "exp")
hist(col.vdp$condition)
}
