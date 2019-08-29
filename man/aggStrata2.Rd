\name{aggrStrata2}
\Rdversion{1.5}
\alias{aggrStrata2}
\title{
Builds the "strata" dataframe containing information on target variables Y's
distributions in the different strata, starting from a frame
}
\description{
This function builds the dataframe "strata" considering as input a given domain in the sampling frame. In case a dataframe "model" is given, the anticipated variance in the strata for each target variable is calculated 
}
\usage{
aggrStrata2(dataset, 
            model=NULL, 
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
\item{model}{
	Dataframe with the parameters of the model(s) to be used to calculate
	anticipated variance.
	Default is null.
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
strata <- aggrStrata2(dataset=frame,
                      model=NULL,
                      vett,
                      dominio=1)
}
}
\keyword{ survey }

