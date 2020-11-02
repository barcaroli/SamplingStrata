\name{procBethel}
\Rdversion{1.5}
\alias{procBethel}
\title{
Procedure to apply Bethel algorithm and select a sample from given strata
}
\description{
This function allows to execute a complete procedure from the Bethel optimal allocation to the selection of a sample, without having to optimize the strata, that are supposed to be given and fixed.
}
\usage{
procBethel(framesamp,
           framecens,
           errors,
           sampling_method=c("srs","systematic","spatial"),
           minnumstrat=2) 
}
\arguments{
  \item{framesamp}{
Dataframe containing sampling frame units.
}
  \item{framecens}{
Dataframe containing frame units that must be selected.
}
  \item{errors}{
Dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation related to target variables of the survey.
}
  \item{sampling_method}{
  Parameter for choosing the selection method: "srs", "systematic" and "spatial".
}
  \item{minnumstrat}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
}

}

\value{
List containing 
(1) the selected sample, 
(2) the strata with allocated sampling units,
(3) the take-all strata
(4) the sampling frame with the labels linking to (2)
(5) the frame with take-all units, with the labels linking to (3) 
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{
data("swissmunicipalities")
swissmun <- swissmunicipalities[swissmunicipalities$REG < 4,
                                c("REG","COM","Nom","HApoly",
                                  "Surfacesbois","Surfacescult",
                                  "Airbat","POPTOT")]
ndom <- length(unique(swissmun$REG))
cv <- as.data.frame(list(DOM=rep("DOM1",ndom),
                         CV1=rep(0.10,ndom),
                         CV2=rep(0.10,ndom),
                         domainvalue=c(1:ndom) ))
cv
swissmun$HApoly.cat <- var.bin(swissmun$HApoly,15)
swissmun$POPTOT.cat <- var.bin(swissmun$POPTOT,15)
frame <- buildFrameDF(df = swissmun,
                         id = "COM",
                         X = c("POPTOT.cat","HApoly.cat"),
                         Y = c("Airbat","Surfacesbois"),
                         domainvalue = "REG")
summary(frame)
#----Selection of units to be censused from the frame
ind_framecens <- which(frame$X1 > 9)
framecens <- frame[ind_framecens,]
#----Selection of units to be sampled from the frame
# (complement to the previous)
framesamp <- frame[-ind_framecens,]

a <- procBethel(framesamp,framecens,errors=cv,sampling_method="srs",minnumstrat=2)
head(a$sample)
expected_CV(a$strata)
}
}
\keyword{ survey }

