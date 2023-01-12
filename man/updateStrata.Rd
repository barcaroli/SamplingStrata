\name{updateStrata}
\Rdversion{1.3}
\alias{updateStrata}
\title{
Assigns new labels to atomic strata on the basis of the optimized aggregated strata 
}
\description{
Once optimal stratification has been obtained ('outstrata'), then we need to attribute 
new strata labels to each atomic stratum.
By executing this function, a new dataframe "newstrata" will be obtained with the same structure of the old, 
("strata") but with the addition of a new stratum label.
By indicating "YES" to "writeFile" parameter, the dataframe "newstrata" will be written to a delimited file 
("newstrata.txt").
Also a second delimited file ("strata_aggregation.txt") will be outputted, 
containing the indication of the relations bewteen atomic and aggregated strata. 
}
\usage{
updateStrata(strata, solution, writeFiles = FALSE, outputFolder = file.path(getwd(),"output"))
}
\arguments{
  \item{strata}{
  This is the (mandatory) dataframe containing the information related to the atomic strata to which the 
  optimisation has been applied to. 
}
  \item{solution}{
  List obtained by the execution of the "optimizeStrata" function. The first element of the list is the vector
  of the indices corresponding to the optimal solution.
}
  \item{writeFiles}{
  Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
  Default is "FALSE".
}
	\item{outputFolder}{
A path to a directory ot be used for the outputs. Default is "output"
	}
}

\value{
A dataframe containing the strata
}

\author{
Giulio Barcaroli
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
# updating sampling strata with new strata labels
newstrata <- updateStrata(swissstrata, solution, writeFiles = TRUE)
}
}
\keyword{ survey }

