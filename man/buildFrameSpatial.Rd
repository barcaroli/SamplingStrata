\name{buildFrameSpatial}
\Rdversion{1.3}
\alias{buildFrameSpatial}
\title{
Builds the "sampling frame" dataframe from a dataset containing information all the units in the population of reference including spatial
}
\description{
This function allows to build the information regarding the sampling frame of the population of reference.
Mandatory variables are:
(i) the name of the dataset containing the sampling frame of the population of reference
(ii) an identifier (Id)
(iii) a set of auxiliary variables 'X'
(iv) a set of target variables 'Y'
(v) a set of prediction errors variables 'variance'
(vi) longitude
(vii) latitude
(viii) the indicator of the domain to which the unit belongs
}
\usage{
buildFrameSpatial(df, id, X, Y, variance, lon, lat, domainvalue)
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
	\item{variance}{
  A character vector containing the names of the prediction error variables
  in the frame dataset
	}
	\item{lon}{
  Longitude of the unit
	}
	\item{lat}{
  Latitude of the unit
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
library(sp)
library(gstat)
library(automap)
library(SamplingStrata)
data("meuse")
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")
#################
# kriging
#################
v <- variogram(lead ~ dist + soil, data=meuse)
fit.vgm <- autofitVariogram(lead ~ elev + soil, meuse, model = "Exp")
plot(v, fit.vgm$var_model)
fit.vgm$var_model
g <- NULL
g <- gstat(g, "Pb", lead ~ dist + soil, meuse)
g
vm <- variogram(g)
vm.fit <- fit.lmc(vm, g, vgm(psill=fit.vgm$var_model$psill[2], 
                             model="Exp", range=fit.vgm$var_model$range[2], 
                             nugget=fit.vgm$var_model$psill[1]))
# Prediction on the whole grid
preds <- predict(vm.fit, meuse.grid)
names(preds)
# [1] "Pb.pred" "Pb.var"
preds$Pb.pred <- ifelse(preds$Pb.pred < 0,0,preds$Pb.pred)
df <- NULL
df$id <- meuse.grid$id
df$Pb.pred <- preds@data$Pb.pred
df$Pb.var <- preds@data$Pb.var
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
df$dom1 <- 1
df <- as.data.frame(df)
frame <- buildFrameSpatial(df=df,
                      id="id",
                      X=c("Pb.pred"),
                      Y=c("Pb.pred"),
                      variance=c("Pb.var"),
                      lon="lon",
                      lat="lat",
                      domainvalue = "dom1")
head(frame)
}
}
\keyword{ survey }


