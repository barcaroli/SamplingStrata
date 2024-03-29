\name{selectSample}
\Rdversion{1.3}
\alias{selectSample}
\title{
Selection of a stratified sample from the frame with srswor method 
}
\description{
Once optimal stratification has been obtained, 
and a new frame has been 
built by assigning to the units of the old one the new strata labels, it is possible to select 
a stratified sample from the frame 
with the simple random sampling without replacement (srswor) method. 
The result of the execution of "selectSample" function is a dataframe containing the selected
units, with their weights (inverse of the probabilities of inclusion).
It is possible to output this dataframe in a .csv file.
One more .csv file is produced ("sampling_check"), containing coeherence checks between 
(a) population in frame strata
(b) population in optimised strata
(c) planned units to be selected in optimised strata
(d) actually selected units
(e) sum of weights in each stratum

}
\usage{
selectSample(frame, outstrata, writeFiles = FALSE,outputFolder = file.path(getwd(),"output"), verbatim=TRUE)
}
\arguments{
  \item{frame}{
  This is the (mandatory) dataframe containing the sampling frame, as it has been modified 
  by the execution of the "updateFrame" function. 
}
  \item{outstrata}{
  This is the (mandatory) dataframe containing the information related to resulting   stratification obtained by the execution of "optimizeStrata" function. It should coincide with 'solution$aggr_strata'.
}
  \item{writeFiles}{
  Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
  Default is "FALSE".
}
	\item{outputFolder}{
A path to a directory ot be used for the outputs. Default is "output"
	}
  \item{verbatim}{
  Indicates if information on the drawn sample must be printed or not.
  Default is "TRUE".
}
}

\value{
A dataframe containing the sample
}

\author{
Giulio Barcaroli with contribution from Diego Zardetto
}

 
\examples{
\dontrun{
library(SamplingStrata)
data(swisserrors)
data(swissstrata)
# optimisation of sampling strata
solution <- optimizeStrata ( 
    errors = swisserrors, 
    strata = swissstrata 
)
# updating sampling strata with new strata labels
newstrata <- updateStrata(swissstrata, solution)
# updating sampling frame with new strata labels
data(swissframe)
framenew <- updateFrame(frame=swissframe,newstrata=newstrata)
# selection of sample
sample <- selectSample(frame=framenew,outstrata=solution$aggr_strata)
head(sample)
}
}
\keyword{ survey }

