library(sp)
# locations (155 observed points)
data("meuse")
# grid of points (3103)
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")
library(automap)
kriging_lead = autoKrige(log(lead) ~ dist, meuse, meuse.grid)
plot(kriging_lead,sp.layout = NULL, justPosition = TRUE)
kriging_zinc = autoKrige(log(zinc) ~ dist, meuse, meuse.grid)
plot(kriging_zinc, sp.layout = list(pts = list("sp.points", meuse)))
r2_lead <- 1 - kriging_lead$sserr/sum((log(meuse$lead)-mean(log(meuse$lead)))^2)
r2_lead
r2_zinc <- 1 - kriging_zinc$sserr/sum((log(meuse$zinc)-mean(log(meuse$zinc)))^2)
r2_zinc
df <- NULL
df$id <- meuse.grid$id
df$lead.pred <- kriging_lead$krige_output@data$var1.pred
df$lead.var <- kriging_lead$krige_output@data$var1.var
df$zinc.pred <- kriging_zinc$krige_output@data$var1.pred
df$zinc.var <- kriging_zinc$krige_output@data$var1.var
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
df$dom1 <- 1
df <- as.data.frame(df)
head(df)
library(SamplingStrata)
# Load the local overrides of the spatial-method functions (must be sourced
# after library(SamplingStrata) so assignInNamespace() below can find the
# loaded namespace to patch).
source("aggrStrata_spatial.r")
source("buildFrameSpatial.R")
source("strataGenalgSpatial.R")
source("optimizeStrataSpatial.R")
source("KmeansSolutionSpatial.R")
source("patchSpatialMethod.R")
frame <- buildFrameSpatial(df=df,
                           id="id",
                           X=c("lead.pred","zinc.pred"),
                           Y=c("lead.pred","zinc.pred"),
                           variance=c("lead.var","zinc.var"),
                           lon="lon",
                           lat="lat",
                           domainvalue = "dom1")

cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.01,1),
                         CV2=rep(0.01,1),
                         domainvalue=c(1:1) ))
set.seed(1234)
km <- KmeansSolutionSpatial(frame,
                            errors = cv,
                            fitting = c(r2_lead,r2_zinc),
                            range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
                            kappa=1,
                            nstrata=NA,
                            maxclusters = 5) 


set.seed(1234)
# system.time(
solution <- SamplingStrata::optimStrata (
  method = "spatial",
  errors=cv,
  framesamp=frame,
  iter = 25,
  pops = 10,
  nStrata = 5,
  fitting = c(r2_lead,r2_zinc),
  range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
  # fitting = 0 is rejected with an explanatory error (see buildStrataDFSpatial.cpp
  # and aggrStrata_spatial.r): it means "the kriging model explains none of the
  # variance", not "no spatial correlation". To disable the spatial-correlation
  # adjustment instead, use fitting = c(1,1) together with range = c(0,0).
  # NOTE: the buildStrataDFSpatial fix requires SamplingStrata to be rebuilt from
  # the fixed src/buildStrataDFSpatial.cpp and reinstalled to take effect here.
  kappa=1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = FALSE)
# )
# detach("package:SamplingStrata", unload = TRUE)
# library(SamplingStrataPlus)
# set.seed(1234)
# system.time(
#   solution2 <- SamplingStrataPlus::optimStrata (
#     method = "spatial",
#     errors=cv,
#     framesamp=frame,
#     iter = 20,
#     pops = 10,
#     nStrata = 5,
#     fitting = c(r2_lead,r2_zinc),
#     range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
#     kappa=1,
#     writeFiles = FALSE,
#     showPlot = TRUE,
#     parallel = FALSE)
# )
framenew <- solution$framenew
outstrata <- solution$aggr_strata
plotStrata2d(framenew,outstrata,domain=1,vars=c("X1","X2"),
             labels=c("Lead","Zinc"))

frameres <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres2 <- SpatialPixelsDataFrame(points=frameres[c("LON","LAT")], data=framenew)
frameres2$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres2,c("LABEL"), col.regions=bpy.colors(5))

s <- selectSampleSpatial(framenew,outstrata,coord_names=c("LON","LAT"))

coordinates(s) <- ~LON+LAT
proj4string(s) <- CRS("+init=epsg:28992")
s$LABEL <- as.factor(s$LABEL)
library(mapview)
mapview(s,zcol="LABEL", map.types = c("OpenStreetMap"))

