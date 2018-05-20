\name{swissstrata}
\Rdversion{1.1}
\alias{swissstrata}
\docType{data}
\title{
Dataframe containing information on strata in the swiss municipalities frame
}
\description{
Dataframe containing information on strata in the swiss municipalities frame
}
\usage{data(swissframe)}
\format{
  The "swissstrata" dataframe contains a row per stratum with the following variables: 
  \describe{
    \item{STRATO}{
	Identifier of the stratum (character)
	}
    \item{N}{
	Number of population units in the stratum (numeric) 
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
    \item{M1}{
	Mean in the stratum of Y1 (number of men and women aged between 0 and 19)(numeric)
	}
    \item{M2}{
	Mean in the stratum of Y2 (number of men and women aged between 20 and 39) (numeric)
	}
    \item{M3}{
	Mean in the stratum of Y3 (number of men and women aged between 40 and 64) (numeric)
	}
	\item{M4}{
	Mean in the stratum of Y4 (number of men and women aged between 64 and over) (numeric)
	}
    \item{S1}{
	Standard deviation in the stratum of Y1 (number of men and women aged between 0 and 19)(numeric)
	}
    \item{S2}{
	Standard deviation in the stratum of Y2 (number of men and women aged between 20 and 39) (numeric)
	}
    \item{S3}{
	Standard deviation in the stratum of Y3 (number of men and women aged between 40 and 64) (numeric)
	}
	\item{S4}{
	Standard deviation in the stratum of Y4 (number of men and women aged between 64 and over) (numeric)
	}
    \item{cens}{
	Flag (1 indicates a take all straum, 0 a sampling stratum) (numeric) 
	Default = 0
	}
    \item{cost}{
	Cost per interview in each stratum.  Default = 1 (numeric)
	}
    \item{DOM1}{
	Value of domain to which the stratum belongs Default = 1 (factor or numeric)
	}
  }
}

\examples{
# data(swissstrata)
# head(swissstrata)
}
\keyword{datasets}
