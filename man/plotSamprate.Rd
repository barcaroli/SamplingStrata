\name{plotSamprate}
\Rdversion{1.3}
\alias{plotSamprate}
\title{
Plotting sampling rates in the different strata for each domain in the solution.
}
\description{
Once the optimization step has been carried out, by applying this function it is
possible to obtain the visualization of the proportion of sampling units in the
different strata for each domain in the obtained solution.
}
\usage{
plotSamprate(solution, dom)
}
\arguments{
  \item{solution}{
Solution obtained by executing optimizeStrata
	}
  \item{dom}{
Identification of the domain
	}
}
\value{
Plot
}
 
\examples{
\dontrun{
library(SamplingStrata)
data(swisserrors)
data(swissstrata)
# optimisation of sampling strata
solution <- optimizeStrata ( 
    errors = swisserrors, 
    strata = swissstrata, 
)
# plot of the sampling rates in strata of domain 1
plotSamprate(solution, 1)
}
}
\keyword{ survey }

