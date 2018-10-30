\name{updateFrame2}
\Rdversion{1.3-1}
\alias{updateFrame2}
\title{
Updates the initial frame on the basis of the optimized stratification (only after 'optimizeStrata2') 
}
\description{
After the execution of 'optimizeStarat2, once optimal stratification has been obtained
it is important to report the new classification of units in the 
sampling frame by attributing new strata labels to each unit.
By executing this function, a new frame will be obtained with the same structure of the old, 
but with the addition of a new stratum label.
 
}
\usage{
updateFrame(frame, solution, writeFiles = FALSE)
}
\arguments{
  \item{frame}{
  This is the (mandatory) dataframe containing the sampling frame. 
}
  \item{solution}{
  The solution obtained by the execution of the 'optimizeStrata2' function. 
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
solution <- optimizeStrata2 ( 
    errors = swisserrors, 
    strata = swissstrata)
# updating sampling frame with new strata labels
framenew <- updateFrame2(frame=swissframe, solution, writeFiles = FALSE)
}
}
\keyword{ survey }

