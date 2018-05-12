\name{updateFrame}
\Rdversion{1.1}
\alias{updateFrame}
\title{
Updates the initial frame on the basis of the optimized stratification 
}
\description{
Once optimal stratification has been obtained, and new labels have been attributed to initial
atomic strata ("newstrata"), it is important to report the new classification of units in the 
sampling frame by attributing new strata labels to each unit.
By executing this function, a new frame will be obtained with the same structure of the old, 
but with the addition of a new stratum label.
The initial frame must contain a variable named 'domainvalue' that indicates the 
same values of the domain that has been used with the 'optimizeStrata' function. 
If no domains have been defined, this variable will contains all 1's, but it must exist
 
}
\usage{
updateFrame(frame, newstrata, writeFiles = FALSE)
}
\arguments{
  \item{frame}{
  This is the (mandatory) dataframe containing the sampling frame. 
}
  \item{newstrata}{
  This is the (mandatory) dataframe containing the information related to the optimisation applied to 
  initial stratification (new labels applied to atomic strata). It is produced by executing 
  the "updateStrata" function. 
}
  \item{writeFiles}{
  Flag to write or not the new sampling frame into the working directory.
  Default is "FALSE"
}
}

\value{
A dataframe containing the frame
}

\author{
Giulio Barcaroli
}

\examples{
#
# The following example is realistic, but is time consuming
#
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
    initialStrata = as.numeric(table(swissstrata$DOM1)),
    addStrataFactor = 0.01,    
    minnumstr = 2,
    iter = 60,
    pops = 20,
    mut_chance = 0.05,
    elitism_rate = 0.2,
    highvalue = 100000000,
    suggestions = NULL,
	realAllocation = TRUE,
	writeFiles = FALSE,
	showPlot = TRUE)
# updating sampling strata with new strata labels
newstrata <- updateStrata(swissstrata, solution, writeFiles = TRUE)
# updating sampling frame with new strata labels
data(swissframe)
framenew <- updateFrame(frame=swissframe, newstrata=newstrata, writeFiles = TRUE)
}
}
\keyword{ survey }

