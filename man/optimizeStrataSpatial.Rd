\name{optimizeStrataSpatial}
\Rdversion{1.3-1}
\alias{optimizeStrataSpatial}
\title{
Best stratification of a sampling frame for multipurpose surveys considering also spatial correlation
}
\description{
This function runs a set of other functions to optimise the stratification of a sampling frame, only when stratification variables are of the continuous type, and if there is also 
a component of spatial autocorrelation in frame units.
}
\usage{
	optimizeStrataSpatial(
            errors, 
            framesamp,
            framecens = NULL, 
            strcens = FALSE, 
            alldomains = TRUE, 
            dom = NULL, 
            nStrata = 5, 
            fitting=1,
            range=1,
            kappa=3,
            minnumstr = 2, 
            iter = 50, 
            pops = 20, 
            mut_chance = NA, 
            elitism_rate = 0.2, 
            highvalue = 1e+08, 
            suggestions = NULL, 
            realAllocation = TRUE, 
            writeFiles = FALSE, 
            showPlot = TRUE, 
            parallel = TRUE, 
            cores
	)
}
\arguments{
  \item{errors}{
  This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
  related to target variables of the survey.
}
  \item{framesamp}{
  This is the (mandatory) dataframe containing the information related to the sampling frame.
}
  \item{framecens}{
  This the (optional) dataframe containing the units to be selected in any case. It has same structure than "frame" dataframe.
}
  \item{strcens}{
  Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
}
  \item{alldomains}{
  Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
  then a value must be given to parameter 'dom'.
}
  \item{dom}{
  Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval 
  (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
}
  \item{nStrata}{
  Indicates the number of strata for each variable.
}
\item{fitting}{
	Fitting of the model(s).
	Default is 1.
}
\item{range}{
	Maximum range for spatial autocorrelation. It is a vector with as many elements as the number of target variables Y.
}	
\item{kappa}{
	Factor used in evaluating spatial autocorrelation. 
	Default is 3.
}	
\item{minnumstr}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
}
  \item{iter}{
  Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
}
  \item{pops}{
  The dimension of each generations in terms of individuals. Default is 20.
}
  \item{mut_chance}{
  Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
  of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
  but a slower convergence, while low values permit a faster convergence, but the final solution can be
  distant from the optimal one. Default is NA, in correspondence of which it is computed as 1/(vars+1) 
  where vars is the length of elements in the solution.
}
  \item{elitism_rate}{
  This  parameter indicates the rate of better solutions that must be preserved from one generation
	to another. Default is 0.2 (20%).
}
  \item{highvalue}{
  Parameter for genetic algorithm. In should not be changed
}
  \item{suggestions}{
  Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL. 
}
  \item{realAllocation}{
  If FALSE, the allocation is based on INTEGER values; if TRUE, the allocation is based on REAL values.
  Default is TRUE.
}  
  \item{writeFiles}{
  Indicates if the various dataframes and plots produced during the execution have to be written in the working
  directory.
  Default is FALSE.
}
  \item{showPlot}{
  Indicates if the plot showing the trend in the value of the objective function has to be shown or not. In parallel = TRUE, this defaults to FALSE
  Default is TRUE.
}
  \item{parallel}{
  Should the analysis be run in parallel.
  Default is TRUE.
}
  \item{cores}{
  If the analysis is run in parallel, how many cores should be used.
  If not specified n-1 of total available cores are used OR if number of domains < (n-1) cores, then number of cores equal to number of domains are used.
}
}

\value{
A list containing 
(1) the vector of the solution, 
(2) the optimal aggregated strata,
(3) the total sampling frame with the label of aggregated strata
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
#################
# Optimization
#################
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
set.seed(1234)
solution <- optimizeStrataSpatial (
  errors=cv, 
  framesamp=frame,
  iter = 50,
  pops = 10,
  nStrata = 3,
  fitting = 1,
  range = fit.vgm$var_model$range[2],
  kappa=1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE)
framenew <- solution$framenew
outstrata <- solution$aggr_strata
############################
# Analysis and visualization
############################
sum(round(outstrata$SOLUZ))
expected_CV(outstrata)
ss <- summaryStrata(framenew,outstrata)
frameres <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres2 <- SpatialPixelsDataFrame(points=frameres[c("LON","LAT")], data=framenew)
frameres2$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres2,c("LABEL"), col.regions=bpy.colors(5))
}
}
\keyword{ survey }

