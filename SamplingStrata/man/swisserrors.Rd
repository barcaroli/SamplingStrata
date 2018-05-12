\name{swisserrors}
\Rdversion{1.1}
\alias{swisserrors}
\docType{data}
\title{
Precision constraints (maximum CVs) as input for Bethel allocation
}
\description{
Dataframe containing precision levels (expressed in terms of acceptable CV's)
}
\usage{data(errors)}
\format{
  The constraint data frame (swisserrors) contains a row per each domain value with the following variables: 
  \describe{
    \item{DOM}{
    Type of domain code (factor)
    } 
    \item{CV1}{
    Planned coefficient of variation for first variable Y1 (number of men and women aged between 0 and 19) (numeric)
    } 
    \item{CV2}{
    Planned coefficient of variation for second variable Y2 (number of men and women aged between 20 and 39) (numeric)
    }
    \item{CV3}{
    Planned coefficient of variation for third variable Y3 (number of men and women aged between 40 and 64) (numeric)
    } 
    \item{CV4}{
    Planned coefficient of variation for forth variable Y4 (number of men and women aged between 65 and over) (numeric)
    } 	
    \item{domainvalue}{
    Value of the domain to which the constraints refer (numeric)
    }
  }
}

\examples{
## data(swisserrors)
## swisserrors
}
\keyword{datasets}
