\name{buildErrorsDF}
\Rdversion{1.5}
\alias{buildErrorsDF}
\title{
Builds the "errors" dataframe containing information on the precision constraints set on the target variables
}
\description{
This function allows to build the dataframe containing the precision constraints to be set on the target variables 
of the sampling survey. It can be applied to the frame dataframe or to the strata dataframe, once they have been 
defined. The dataframe values of the precision constraints will not be initialized, so the user has to do it.
}
\usage{
buildErrorsDF(frame=NULL,
              strata=NULL)
}
\arguments{
  \item{frame}{
  This is the name of the dataframe containing the sample data, or the frame data. 
  It is strictly required that auxiliary information is organised in variables named 
  as X1, X2, ... , Xm (there should be at least one of them) and the target variables 
  are denoted by Y1, Y2, ... , Yn.
  In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe, 
  containing   the weigths associated to each sampling unit
	}
\item{strata}{
  In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
  to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
  target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters 
  beta, sig2, gamma (> 1 in case of heteroscedasticity). 
  Default is NULL.
}
}

\value{
A dataframe containing precision constraints on the target variables.
}
\author{
Giulio Barcaroli
}

\examples{
\dontrun{
library(SamplingStrata)
data(swissmunicipalities)
swissmunicipalities$HApoly.cat <- var.bin(swissmunicipalities$HApoly,15)
swissmunicipalities$POPTOT.cat <- var.bin(swissmunicipalities$POPTOT,15)
frame <- buildFrameDF(df = swissmunicipalities,
                       id = "COM",
                       X = c("POPTOT.cat","HApoly.cat"),
                       Y = c("Airbat","Surfacesbois"),
                       domainvalue = "REG")
strata <- buildStrataDF(frame)
ndom <- length(unique(swissmunicipalities$REG))
cv <- buildErrorsDF(strata=strata)
#cv <- buildErrorsDF(frame=frame)
cv$CV1 <- c(rep(0.03,7))
cv$CV2 <- c(rep(0.05,7))
}
}
\keyword{ survey }

