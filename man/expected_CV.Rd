\name{expected_CV}
\Rdversion{1.3-1}
\alias{expected_CV}
\title{
Expected coefficients of variation of target variables Y
}
\description{
Once optimized the sampling frame, this function allows to calculate the expected coefficients of variation on the aggregated strata of the current optimized solution, and to compare them to the precision constraints.
}
\usage{
expected_CV(strata)
}
\arguments{
  \item{strata}{
Aggregated strata in the solution obtained by the execution of the 'optimizedStrata' function
	}
}
\value{
Matrix containing values of the CV's in the different domains
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
# calculate CV's on Y's
expected_CV(solution$aggr_strata)
# compare to precision constraints
swisserrors
}
}
\keyword{ survey }

