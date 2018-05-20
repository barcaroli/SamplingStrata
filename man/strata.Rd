\name{strata}
\Rdversion{1.1}
\alias{strata}
\docType{data}
\title{
Dataframe containing information on strata in the frame
}
\description{
Dataframe containing information on strata in the frame
}
\usage{data(strata)}
\format{
  The strata data frame (strata) contains a row per stratum with the following variables: 
  \describe{
    \item{stratum}{
	Identifier of the stratum (numeric)
	}
    \item{N}{
	Number of population units in the stratum (numeric) 
	}
    \item{X1}{
	Value of first auxiliary variable X1 in the stratum (factor)
	}
    \item{Xi}{
	Value of i-th auxiliary variable Xi in the stratum (factor)
	}
    \item{Xk}{
	Value of last auxiliary variable Xk in the stratum (factor)
	}
    \item{M1}{
	Mean in the stratum of the first variable Y1 (numeric)
	}
    \item{Mj}{
	Mean in the stratum of the j-th variable Yt (numeric)
	}
    \item{Mn}{
	Mean in the stratum of the last variable Y (numeric)
	}
    \item{S1}{
	Standard deviation in the stratum of the first variable Y (numeric)
	}
    \item{Sj}{
	Standard deviation in the stratum of the j-th variable Yt (numeric)
	}
    \item{Sn}{
	Standard deviation in the stratum of the last variable Y (numeric)
	}
    \item{cens}{
	Flag (1 indicates a take all straum, 0 a sampling stratum) (numeric) 
	Default = 0
	}
    \item{cost}{
	Cost per interview in each stratum.  Default = 1 (numeric)
	}
    \item{DOM1}{
	Value of domain to which the stratum belongs (factor or numeric)
	}
  }
}
\details{
Note: the names of the variables must be the ones indicated above
}

\examples{
# data(strata)
# head(strata)
}
\keyword{datasets}
