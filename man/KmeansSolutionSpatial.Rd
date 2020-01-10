\name{KmeansSolutionSpatial}
\Rdversion{1.3}
\alias{KmeansSolutionSpatial}
\title{
Initial solution obtained by applying kmeans clustering of frame units
}
\description{
This function has to be used only in conjunction with "optimizeStrataSpatial", i.e. in the
case of optimizing with only continuous stratification variables
and with a component of spatial autocorrelation.
The function "KmeansSolutionSpatial" has a twofold objective:
- to give indications about a possible best number of final strata (by fixing a 
convenient value for "maxclusters", and leaving NA to "nstrata";
- to give an initial solution fo the optimization step.
If the parameter "nstrata" is not indicated, the optimal number of clusters 
is determined inside each domain, and the overall solution is obtained by concatenating optimal 
clusters obtained in domains.
The result is a dataframe with two columns: the first indicates the clusters, 
the second the domains.
}
\usage{
KmeansSolutionSpatial(frame, 
                    fitting = 1, 
                    range = c(0), 
                    kappa = 3, 
                    errors, 
                    nstrata = NA, 
                    minnumstrat = 2, 
                    maxclusters = NA, 
                    showPlot = TRUE)
}
\arguments{
  \item{frame}{
  The (mandatory) dataframe containing the information related to each unit in the sampling frame.
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
\item{errors}{
  The (mandatory) dataframe containing the precision constraints on target variables.
}
\item{nstrata}{
  Number of aggregate strata (if NULL, it is optimized by varying the number of cluster
  from 2 to half number of atomic strata). Default is NA.
}
\item{minnumstrat}{
  Minimum number of units to be selected in each stratum. Default is 2.
}
\item{maxclusters}{
Maximum number of clusters to be considered in the 
execution of kmeans algorithm. If not indicated it will be set equal to the
number of atomic strata divided by 2.
}
\item{showPlot}{
Allows to visualise on a plot the different sample sizes for each number of aggregate strata. Default is TRUE.
}
}

\value{
A dataframe containing the solution
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
#################
# data
#################
# locations (155 observed points)
data("meuse")
# grid of points (3103)
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
# model    psill    range
# 1   Nug 1524.895   0.0000
# 2   Exp 8275.431 458.3303
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
df$Pb.pred <- preds@data$Pb.pred
df$Pb.var <- preds@data$Pb.var
df$dom1 <- 1
df <- as.data.frame(df)
df$id <- meuse.grid$id
#####################################
# Optimization with kmeans clustering
#####################################
frame <- buildFrameDF(df=df,
                      id="id",
                      X=c("Pb.pred"),
                      Y=c("Pb.pred"),
                      domainvalue = "dom1")
frame$var1 <- df$Pb.var
frame$lon <- meuse.grid$x
frame$lat <- meuse.grid$y
cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.05,1),
                         domainvalue=c(1:1) ))
km <- KmeansSolutionSpatial(frame,
                            errors = cv,
                            fitting = 1,
                            range = fit.vgm$var_model$range[2],
                            kappa=1,
                            nstrata=NA,
                            maxclusters = 5)
############################
# Analysis and visualization
############################
strataKm <- aggrStrataSpatial(dataset=frame,
                        fitting = 1,
                        range = fit.vgm$var_model$range[2],
                        kappa=1,
                        vett=km$suggestions,
                        dominio=1)
strataKm$SOLUZ <- bethel(strataKm,cv)
sum(strataKm$SOLUZ)
framenew <- frame
framenew$LABEL <- km$suggestions
strataKm$STRATO <- strataKm$stratum
ssKm <- summaryStrata(framenew,strataKm)
ssKm
frameres <- SpatialPointsDataFrame(data=framenew, 
                                   coords=cbind(framenew$lon,framenew$lat) )
frameres2 <- SpatialPixelsDataFrame(points=frameres[c("lon","lat")], 
                                    data=framenew)
frameres2$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres2,c("LABEL"), col.regions=bpy.colors(5))
}
}
\keyword{ survey }

