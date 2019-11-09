\name{computeGamma}
\Rdversion{1.5}
\alias{computeGamma}
\title{
Function that allows to calculate a heteroscedasticity index,
together with associate prediction variance,
to be used by the optimization step 
to correctly evaluate the standard deviation in the strata
due to prediction errors.
}
\description{
When the anticipated variance has to be calculated during the execution
of the optimization step, his function allows to calculate a heteroscedasticity index,
together with associate prediction variance,
to be used to correctly evaluate the variance in the strata.
The function returns a list where the first object is the heteroscedasticity index and the 
second is the associated standard deviation in the strata
due to prediction errors. 
The two parameters are calculated in this way:
(i) residuals 'e' are grouped in clusters defined by values of the explanatory variable 'x';
(ii) a model is fitted by considering log(e) and log(mean(x)) values;
(iii) the intercept is the value of standard deviation of residuals;
(iv) the slope is the value of the heteroscedasticity index.
These two values can be passed as parameters of the model, or used to calculate
prediction errors for ach unit in the frame.
}
\usage{
computeGamma(e,x,nbins,showPlot)
}
\arguments{
  \item{e}{
  This is the variable that contains prediction errors (residuals) of the model.
	}
\item{x}{
  This is the variable that contains explanatory variable in the model.
}
\item{nbins}{
  Number of bins to be passed to the 'var.bin' function. Default is 6.
}
\item{showPlot}{
  Visualization of plots. Default is TRUE.
}
}
\value{
A list containing: (i) the value of the heteroscedasticity index, (ii) associated standard deviation, (iii) R^2 of the interpolating model.
}
\author{
Marco Ballin, Giulio Barcaroli 
}

\examples{
\dontrun{
data("swissmunicipalities")
swiss_sample <- swissmunicipalities[sample(c(1:nrow(swissmunicipalities)),500),]
mod_Airind_POPTOT <- lm(swiss_sample$Airind ~ swiss_sample$POPTOT)
computeGamma(mod_Airind_POPTOT$residuals,
             swiss_sample$POPTOT,
             nbins = 8)
#     gamma     sigma  r.square 
# 0.8029292 0.0150446 0.9598539
}
}
\keyword{ survey }

