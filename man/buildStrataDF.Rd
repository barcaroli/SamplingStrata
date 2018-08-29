\name{buildStrataDF}
\Rdversion{1.3}
\alias{buildStrataDF}
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
buildStrataDF(dataset,model=NULL,progress=TRUE)
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
  In case the Y variables are not directly observed, but are estimated by means of other explicative variables, in order
  to compute the anticipated variance, information on models are given by a dataframe "model" with as many rows as the
  target variables. Each row contains the indication if the model is linear o loglinear, and the values of the model parameters 
  beta, sig2, gamma (> 1 in case of heteroscedasticity). 
  Default is NULL.
}
\item{progress}{
  If set to TRUE, a progress bar is visualised during the execution. Default is TRUE.
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
# Plain example without model
data(swissframe)
strata <- buildStrataDF(dataset=swissframe,model=NULL)
head(strata)

# More complex example with models
library(SamplingStrata)
data(swissmunicipalities)
swiss <- swissmunicipalities[,c("HApoly","Surfacesbois","Airind","POPTOT")]
Y1 = swiss$Surfacesbois
X1 = swiss$HApoly
mod1 <- lm( Y1 ~ X1 )
summary(mod1)
mod1$coefficients[2]
summary(mod1)$sigma

Y2 = swiss$Airind
X2 = swiss$POPTOT
plot(log(X2[X2>0]),log(Y2[X2>0]))
mod2 <- lm( log(Y2[X2 > 0 & Y2>0]) ~ log(X2[X2 > 0 & Y2>0]) )
summary(mod2)
mod2$coefficients[2]
summary(mod2)$sigma

swiss$id <- c(1:nrow(swiss))
swiss$dom <- 1
frame <- buildFrameDF(swiss,id="id",X="id",Y=c("HApoly","POPTOT"),domainvalue="dom")

model <- NULL
model$type[1] <- "linear"
model$beta[1] <- mod1$coefficients[2]
model$sig2[1] <- summary(mod1)$sigma
model$gamma[1] = 2
model$type[2] <- "loglinear"
model$beta[2] <- mod2$coefficients[2]
model$sig2[2] <- summary(mod2)$sigma
model$gamma[2] = NA
model <- as.data.frame(model)

strata <- buildStrataDF(dataset=frame, model=model)
}
}
\keyword{ survey }

