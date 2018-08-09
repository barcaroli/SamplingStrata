\name{selectSampleSystematic}
\Rdversion{1.3}
\alias{selectSampleSystematic}
\title{
Selection of a stratified sample from the frame with systematic method 
}
\description{
Once optimal stratification has been obtained, 
and a new frame has been 
built by assigning to the units of the old one the new stratum labels (by means of 
"updateFrame" function), it is possible to select a stratified sample from the frame 
with the systematic method, that is a selection that begins selecting the first 
unit by an initial ramndomly chosen starting point, and proceeding in selecting other units 
by adding an interval that is the inverse of the sampling rate in the stratum.

This selection method can be useful if associated to a particular ordering of the selection
frame, where the ordering variable(s) can be considered as additional 
stratum variable(s). 

The result of the execution of "selectSampleSystematic" function is a dataframe containing selected
units, with the probabilities of inclusion.
It is possible to output this dataframe in a .csv file.
One more .csv file is produced ("sampling check"), containing coeherence checks between 
(a) population in frame strata
(b) population in optimised strata
(c) planned units to be selected in optimised strata
(d) actually selected units
(e) sum of weights in each stratum

}
\usage{
selectSampleSystematic(frame, outstrata, sortvariable = NULL, writeFiles = FALSE,verbatim=TRUE)
}
\arguments{
  \item{frame}{
  This is the (mandatory) dataframe containing the sampling frame, as it has been modified 
  by the execution of the "updateFrame" function. 
  Name of stratum variable must be 'strato'.
}
  \item{outstrata}{
  This is the (mandatory) dataframe containing the information related to resulting 
  stratification obtained by the execution of "optimizeStrata" function.
  Name of stratum variable must be 'strato'.
}
  \item{sortvariable}{
  This is the name of the variable to be used as ordering variable inside each stratum before proceeding to the systematic selection. It must be previously added to the selection frame. Default is NULL.
}
  \item{writeFiles}{
  Indicates if at the end of the processing the resulting strata will be outputted in a delimited file.
  Default is "FALSE".
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
#
# The following example is realistic, but is time consuming
#
\dontrun{
library(SamplingStrata)
data(swisserrors)
data(swissstrata)
# optimisation of sampling strata
solution <- optimizeStrata ( 
    errors = swisserrors, 
    strata = swissstrata)
# updating sampling strata with new strata labels
newstrata <- updateStrata(swissstrata, solution)
# updating sampling frame with new strata labels
data(swissframe)
framenew <- updateFrame(frame=swissframe,newstrata=newstrata)
# adding variable "POPTOT" to framenew
data("swissmunicipalities")
framenew <- merge(framenew,swissmunicipalities[,c("REG","Nom","POPTOT")],
                  by.x=c("REG","ID"),by.y=c("REG","Nom"))
# selection of sample with systematic method
sample <- selectSampleSystematic(frame=framenew,
								outstrata=solution$aggr_strata,
								sortvariable="POPTOT")
head(sample)
}
}
\keyword{ survey }

