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
data("swissmunicipalities")
swissmunicipalities = swissmunicipalities[swissmunicipalities$REG==1,]
data("errors")
errors$CV1 <- 0.1
errors$CV2 <- 0.1
errors <- errors[rep(row.names(errors),7),]
errors$domainvalue <- c(1:7)
errors
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
swissmunicipalities$domain = 1
frame <- buildFrameDF(swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("Surfacesbois","Surfacescult"),
                      Y = c("Pop020", "Pop2040")
)
solution <- optimizeStrata2 (
            errors, 
            frame,
            cens = NULL, 
            strcens = FALSE, 
            alldomains = TRUE, 
            dom = NULL, 
            nStrata = 5, 
            minnumstr = 2, 
            iter = 50, 
            pops = 20, 
            mut_chance = NA, 
            elitism_rate = 0.2, 
            highvalue = 1e+08, 
            suggestions = NULL, 
            realAllocation = TRUE, 
            writeFiles = FALSE, 
            showPlot = TRUE, 
            parallel = FALSE) 
framenew <- updateFrame2(frame,solution)
}
}
\keyword{ survey }

