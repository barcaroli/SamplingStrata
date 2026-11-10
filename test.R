library(SamplingStrata)
############################
# Example of "atomic" method
############################
data(swissmunicipalities)
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
frame <- buildFrameDF(df = swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("POPTOT","HApoly"),
                      Y = c("Surfacesbois", "Airind"))
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
                         CV1 = rep(0.1,ndom),
                         CV2 = rep(0.1,ndom),
                         domainvalue = c(1:ndom)))
strata <- buildStrataDF(frame)
kmean <- KmeansSolution(strata,cv,maxclusters=30)
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
# with the atomic method the KmeansSolution() output is used directly
# as suggestions (prepareSuggestion() is only for continuous / spatial)
solution <- optimStrata(method ="atomic",
                        framesamp = frame,
                        errors = cv,
                        nStrata = nstrat,
                        suggestions = kmean,
                        iter = 50,
                        pops = 10)
outstrata <- solution$aggr_strata
framenew <- solution$framenew
s <- selectSample(framenew, outstrata)
################################
# Example of "continuous" method
################################
kmean <- KmeansSolution2(frame = frame,
                         errors = cv,
                         maxclusters = 10)
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
sugg <- prepareSuggestion(kmean = kmean,
                          frame = frame,
                          nStrata = nstrat,
                          seeds = c("kmean", "univariate"))
attr(sugg, "seeds")   
evalSuggestion(frame, errors = cv, sugg, which = 1)
evalSuggestion(frame, errors = cv, sugg, which = 2)
solution <- optimStrata(method = "continuous",
                        framesamp = frame,
                        errors = cv,
                        nStrata = nstrat,
                        iter = 50,
                        pops = 10,
                        suggestions = sugg)
framenew <- solution$framenew
outstrata <- solution$aggr_strata
s <- selectSample(framenew,outstrata)
#############################
# Example of "spatial" method
#############################
library(sp)
data("meuse")
data("meuse.grid")
meuse.grid$id <- c(1:nrow(meuse.grid))
coordinates(meuse) <- c('x','y')
coordinates(meuse.grid) <- c('x','y')
library(gstat)
library(automap)
v <- variogram(lead ~ dist + soil, data = meuse)
fit.vgm.lead <- autofitVariogram(lead ~ dist + soil,
                                 meuse,
                                 model = "Exp")
plot(v, fit.vgm.lead$var_model)
lead.kr <- krige(lead ~ dist + soil, meuse, meuse.grid,
                 model = fit.vgm.lead$var_model)
lead.pred <- ifelse(lead.kr[1]$var1.pred < 0,0, lead.kr[1]$var1.pred)
lead.var <- ifelse(lead.kr[2]$var1.var < 0, 0, lead.kr[2]$var1.var)
df <- as.data.frame(list(dom = rep(1,nrow(meuse.grid)),
                         lead.pred = lead.pred,
                         lead.var = lead.var,
                         lon = meuse.grid$x,
                         lat = meuse.grid$y,
                         id = c(1:nrow(meuse.grid))))
frame <-buildFrameSpatial(df = df,
                          id = "id",
                          X = c("lead.pred"),
                          Y = c("lead.pred"),
                          variance = c("lead.var"),
                          lon = "lon",
                          lat = "lat",
                          domainvalue = "dom")
cv <- as.data.frame(list(DOM = rep("DOM1",1),
                         CV1 = rep(0.05,1),
                         domainvalue = c(1:1) ))
solution <- optimStrata(method = "spatial",
                        errors = cv,
                        framesamp = frame,
                        iter = 25,
                        pops = 10,
                        nStrata = 5,
                        fitting = 1,
                        kappa = 1,
                        range = fit.vgm.lead$var_model$range[2])
framenew <- solution$framenew
outstrata <- solution$aggr_strata
frameres <- SpatialPixelsDataFrame(points = framenew[c("LON","LAT")],
                                   data = framenew)
frameres$LABEL <- as.factor(frameres$LABEL)
spplot(frameres,c("LABEL"), col.regions=bpy.colors(5))
s <- selectSample(framenew,outstrata)