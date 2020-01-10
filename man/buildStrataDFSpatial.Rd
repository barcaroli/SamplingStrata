\name{buildStrataDFSpatial}
\Rdversion{1.3}
\alias{buildStrataDFSpatial}
\title{
Builds the "strata" dataframe containing information on target variables Y's
distributions in the different strata, starting from sample data or from a frame
}
\description{
This function allows to build the information regarding strata in the population required as
an input by the algorithm of Bethel for the optimal allocation.
In order to estimate means and standard deviations for target variables Y's, we need data coming from:
(1) a previous round of the survey whose sample we want to plan;
(2) sample data from a survey with variables that are proxy to the ones we are interested to;
(3) a frame containing values of Y's variables (or proxy variables) for all the population.
In all cases, each unit in the dataset must contain auxiliary information (X's variables)
and also target variables Y's (or proxy variables) values: under these conditions it is possible 
to build the dataframe "strata", containing information on the distribution of Y's in the different strata 
(namely, means and standard deviations), together with information on strata (total population, 
if it is to be censused or not, the cost per single interview).
If the information is contained in a sample dataset, a variable named WEIGHT is expected to be 
present. In case of a frame, no such variable is given, and the function will define a WEIGHT variable
for each unit, whose value is always '1'.
Missing values for each Y variable will not be taken into account in the computation of means and standard 
deviations (in any case, NA's can be present in the dataset).
The dataframe "strata" is written to an external file (tab delimited, extension "txt"), and will be
used as an input by the function "optimizeStrata".
}
\usage{
buildStrataDFSpatial(dataset,
                     fitting=c(1),
                     range=c(0),
                     kappa=3,
                     progress=FALSE,
                     verbose=FALSE)
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
\item{progress}{
  If set to TRUE, a progress bar is visualised during the execution. Default is FALSE.
}
\item{verbose}{
  If set to TRUE, information is given about the number of strata generated. Default is FALSE.
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
strata <- buildStrataDFSpatial(dataset=frame,range=800)
}
}
\keyword{ survey }

