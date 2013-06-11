\name{columbus}
\alias{columbus}
\docType{data}
\title{
Columbus crime
}
\description{
Crime rate in planning neighborhoods in Columbus, Ohio in 1980
}
\usage{data(columbus)}
\format{
  A data frame with 49 observations on the following 6 variables.
  \describe{
    \item{\code{houseval}}{a numeric vector}
    \item{\code{income}}{a numeric vector}
    \item{\code{crime}}{a numeric vector}
    \item{\code{distcbd}}{a numeric vector}
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
The data consist of variables for mean housing value, mean household income, residential and vehicle thefts 
combined per thousand people for 1980, distance to the central business district (CBD), and x and y spatial 
coordinates of neighborhood centroids.
}
\source{
Anselin L (1988) Spatial Econometrics: Methods and Models. Kluwer, Dordrecht
}
\references{
Wheeler DC (2007) Diagnostic tools and a remedial method for collinearity in geographically weighted regression.
Environment and Planning A, 39: 2464-2481 
}
\examples{
data(columbus)
plot(columbus$x, columbus$y)
}
\keyword{datasets}
