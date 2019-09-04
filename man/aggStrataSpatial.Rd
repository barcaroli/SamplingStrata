\name{aggrStrataSpatial}
\Rdversion{1.5}
\alias{aggrStrataSpatial}
\title{
Builds the "strata" dataframe containing information on target variables Y's
distributions in the different strata, starting from a frame where units are spatially correlated.
}
\description{
This function builds the dataframe "strata" considering as input a given domain in the sampling frame. The variance in each stratum is calculated by considering also the component of spatial autocorrelation.
}
\usage{
aggrStrataSpatial(dataset,
            fitting,
            range,
            kappa, 
            vett, 
            dominio)
}
\arguments{
  \item{dataset}{
  This is the name of the dataframe containing the sample data, or the frame data. 
  It is strictly required that auxiliary information is organised in variables named 
  as X1, X2, ... , Xm (there should be at least one of them) and the target variables 
  are denoted by Y1, Y2, ... , Yn.
  In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe, 
  containing   the weigths associated to each sampling unit
	}
\item{fitting}{
	Fitting of the model(s).
	Default is 1.
}
\item{range}{
	Maximum range for spatial autocorrelation
}	
\item{kappa}{
	Factor used in evaluating spatial autocorrelation. 
	Default is 3.
}	
\item{vett}{
	vector of values indicating how the units in the dataset must 
	be aggregated in strata.
}	
\item{dominio}{
  Value indicating the domain in the dataset to be processed.
}
}
\value{
A dataframe containing strata
}
\author{
Giulio Barcaroli
}

\examples{
\dontrun{
strata <- aggrStrataSpatial(dataset=frame,
                            fitting=1,
                            range=800,
                            kappa=1)
}
}
\keyword{ survey }

