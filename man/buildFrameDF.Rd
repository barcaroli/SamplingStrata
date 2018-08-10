\name{buildFrameDF}
\Rdversion{1.3}
\alias{buildFrameDF}
\title{
Builds the "sampling frame" dataframe from a dataset containing information on
all the units in the population of reference
}
\description{
This function allows to build the information regarding the sampling frame of the population of reference.
Mandatory variables are:
(i) the name of the dataset containing the sampling frame of the population of reference
(ii) an identifier (Id)
(iii) a set of auxiliary variables X
(iv) a set of target variables Y
(v) the indicator of the domain to which the unit belongs
}
\usage{
buildFrameDF(df, id, X, Y, domainvalue)
}
\arguments{
  \item{df}{
  This is the name of the dataframe containing the information
  on all the units in population of reference. 
	}
  \item{id}{
  This is the name of the identifier in the sampling frame. 
	}
  \item{X}{
  A character vector containing the names of the auxiliary variables
  in the frame dataset
	}
	\item{Y}{
  A character vector containing the names of the target variables
  in the frame dataset
	}
  \item{domainvalue}{
  The name of the variable in the frame dataset that contains the 
  indication of the domains to which the units belong.
	}
}
\value{
A dataframe
}
\author{
Giulio Barcaroli
}

\examples{
\dontrun{
data(swissmunicipalities)
id = "Nom"
X = c("Surfacesbois","Surfacescult")
Y = c("Pop020","Pop2040")
domainvalue = "REG"
frame <- buildFrameDF(swissmunicipalities,id,X,Y,domainvalue)
head(frame)
}
}
\keyword{ survey }


