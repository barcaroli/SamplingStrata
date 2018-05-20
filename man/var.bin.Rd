\name{var.bin}
\Rdversion{1.1}
\alias{var.bin}
\title{
Allows to transform a continuous variable into a categorical ordinal one
by applying a modified version of the k-means clustering function in the 'stats' package.
}
\description{
The optimization of a frame stratification is applicable only in presence of all categorical
auxiliary variables in the frame. If one or more continuous auxiliary variables are in the frame, 
it is necessary to pre-process in order to convert them into categorical (ordinal) variables.
The applied method is the "k-means" clustering method contained in the in "stats" package.
This function ensures that the final result is in an ordered categorical variable.
}
\usage{
var.bin(x,
		bins=3,
		iter.max=100)
}
\arguments{
  \item{x}{
Continuous variable to be transformed into a categorical one 
	}
  \item{bins}{
Number of values of the resulting categorical variable
	}
  \item{iter.max}{
Maximum number of iterations of the clustering algorithm
	}
}
\value{
Binned variable
}
 
\examples{
\donttest{
library(SamplingStrata)
data(swissmunicipalities)
data(swissframe)
swissframe$X1 <- var.bin(swissmunicipalities$POPTOT,bins = 18)
table(swissframe$X1)
tapply(swissmunicipalities$POPTOT,swissframe$X1,mean)
}
}
\keyword{ survey }

