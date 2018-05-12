\name{errors}
\Rdversion{1.1}
\alias{errors}
\docType{data}
\title{
Precision constraints (maximum CVs) as input for Bethel allocation
}
\description{
Dataframe containing precision levels (expressed in terms of acceptable CV's)
}
\usage{data(errors)}
\format{
  The constraint data frame (errors) contains a row per each domain value with the following variables: 
  \describe{
    \item{DOM}{
    Type of domain code (factor)
    } 
    \item{CV1}{
    Planned coefficient of variation for first variable Y1 (numeric)
    } 
    \item{CVj}{
    Planned coefficient of variation for j-th variable Yj (numeric)
    }
    \item{CVn}{
    Planned coefficient of variation for last variable Yn (numeric)
    } 
    \item{domainvalue}{
    Value of the domain to which the constraints refer (numeric)
    }
  }
}
\details{
Note: the names of the variables must be the ones indicated above
}

\examples{
## data(errors)
## errors
}
\keyword{datasets}
