\name{plotSamprate}
\Rdversion{1.1}
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
    cens = NULL,
    strcens = FALSE,
    initialStrata = 3000,
    addStrataFactor = 0.01,    
    minnumstr = 2,
    iter = 60,
    pops = 20,
    mut_chance = 0.05,
    elitism_rate = 0.2,
    highvalue = 100000000,
    suggestions = NULL,
	writeFiles = FALSE)
# plot of the sampling rates in strata of domain 1
plotSamprate(solution, 1)
}
}
\keyword{ survey }

