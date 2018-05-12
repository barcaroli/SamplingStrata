\name{swissframe}
\Rdversion{1.1}
\alias{swissframe}
\docType{data}
\title{
Dataframe containing information on all units in the population of reference
that can be considered as the final sampling unit
(this example is related to Swiss municipalities)
}
\description{
Dataframe containing information on all municipalities in Swiss
(it is a derivation of dataframe "swissmunicipalities" in "sampling" package)
}
\usage{data(swissframe)}
\format{
  The "swissframe" dataframe contains a row per each Swiss municipality with the following variables: 
  \describe{
    \item{progr}{
	Progressive associated to the frame unit (numeric)
	}
	\item{id}{
	Name of the frame unit (character)
	}
    \item{X1}{
	Classes of total population in the municipality (factor with 18 values)
	}
    \item{X2}{
	Classes of wood area in the municipality (factor with 3 values)
	}
    \item{X3}{
	Classes of area under cultivation in the municipality (factor with 3 values)
	}
	\item{X4}{
	Classes of mountain pasture area in the municipality (factor with 3 values)
	}
	\item{X5}{
	Classes of area with buildings in the municipality (factor with 3 values)
	}
	\item{X6}{
	Classes of industrial area in the municipality (factor with 3 values)
	}
    \item{Y1}{
	Number of men and women aged between 0 and 19 (numeric)
	}
    \item{Y2}{
	Number of men and women aged between 20 and 39 (numeric)
	}
    \item{Y3}{
	Number of men and women aged between 40 and 64 (numeric)
	}
	\item{Y4}{
	Number of men and women aged between 65 and over (numeric)
	}
    \item{domainvalue}{
	Value of domain to which the municipality belongs (factor or numeric)
	}
  }
}

\examples{
#data(swissframe)
#head(strata)
}
\keyword{datasets}
