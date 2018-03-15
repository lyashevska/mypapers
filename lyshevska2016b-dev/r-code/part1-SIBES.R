# sequential conditional gaussian simulation on 500 m grid
# model silt+silt2+depth
# binomial and poisson 
# all variables are scaled 

set.seed(5021)

# load libraries
library(gstat)
library(geoR)
library(raster)
library(rgdal)
library(sp)
library(geoRglm)
library(ggplot2)
library(gstat)
library(pscl)
library(MASS)

# load custom functions
source("myfun/myfun.sep.R")
# if nsim=0 then kriging interpolation is done
# back transforming predictions
# christensen version poisson https://github.com/cran/geoRglm/blob/master/R/extensions.R
transf.predP<-function(var1.pred, var1.var){exp(var1.pred+0.5*var1.var)}
transf.varP<-function(var1.pred, var1.var){(exp(2*var1.pred+var1.var))*expm1(var1.var)}
# christensen version binomial https://github.com/cran/geoRglm/blob/master/R/binom.pred.R
### using second order taylor-expansion + facts for N(0,1) [third moment = 0 ; fourth moment = 12].
#plogis(x) == (1+ tanh(x/2))/2
#only for var1.pred<700
transf.predB<-function(var1.pred, var1.var){plogis(var1.pred) + 0.5*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)*var1.var}
transf.varB<-function(var1.pred, var1.var){
    (exp(var1.pred)/(1+exp(var1.pred))^2)^2*var1.var+(1/2)*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)^2*var1.var^2}
# if nsim=1 (or output of glsm.mcmc) then conditional sequential gaussian simulation is performed
# back transforming is done using
antilogit<-function(x){exp(x)/(1+exp(x))} #for bin
#exp for poisson

myscale<-function(x){(x-mean(x))/sd(x)}

##----dat.sep
# read dat.rds
dat<-read.csv("input/dat.csv", header=T)
# remove 4 values in macoma >100
dat<-subset(dat, macoma<100)

dat<- subset(dat, select=c(macoma, silt, depth, x, y))
dat$silt2<-dat$silt^2
dat<-dat[c("macoma", "x", "y","silt", "silt2", "depth")]

# save mean and sd which I use for scaling
meancov<-apply(dat, 2, mean)
sdcov<-apply(dat, 2, sd)
dat<-cbind(dat[,c("macoma", "x", "y")],apply(dat[,c("silt", "silt2", "depth")], 2, myscale))

#datsc<-dat
#write.csv(datsc, "input/datsc.csv")

# model selection
#skip, use ecological model
datP<-dat
datB<-dat
datB$macoma<-ifelse(datB$macoma>0,1,0)

glmP <- glm(formula = macoma~silt+silt2+depth,
	    family = poisson(link=log),
	    data = datP)
summary(glmP)
glmP <- glm(formula = macoma~silt+silt2+depth,
	    family = poisson(link=log),
	    data = datP)
summary(glmP)
glmB <- glm(formula = macoma~silt+silt2+depth,
	    family = binomial,
	    data = datB)
summary(glmB)

#use models to partition data
trendP<-glmP$formula
#macoma ~ mgs + mgs2 + depth
trendB<-glmB$formula
#macoma ~ mgs2 + depth + depth2

dat.sep<-myfun.sep(dat)
table(dat.sep$bin$macoma==0)
table(dat.sep$pois$macoma==0)

# estimate hyperparameters using methods of moment

##----modelvarB
bin<-dat.sep$bin
bin.glm<-glm(formula=trendB, family = "binomial", data = bin)
bin$residB <- resid(bin.glm, type='response')
bin.glm.beta<-as.numeric(bin.glm$coefficients)
coordinates(bin)<-c("x","y")
samplevarB <- variogram(residB ~ 1, data = bin, cutoff=5000)

# fit by eye!
pdf("fig/samplevarB.pdf")
plot(samplevarB)
dev.off()

#
modelvarB <- vgm(psill = 0.02, model = "Sph", range = 2000, nugget = 0.18)
modelvarB <- fit.variogram(samplevarB, model = modelvarB)
linevarB<-variogramLine(modelvarB,maxdist=max(samplevarB$dist))

##----modelvarP
pois<-dat.sep$pois
pois.glm<-glm(formula=trendP, family = "poisson", data = pois)
pois$residP <- resid(pois.glm,type='deviance')
pois.glm.beta<-as.numeric(pois.glm$coefficients)
coordinates(pois)<-c("x","y")
samplevarP <- variogram(residP ~ 1, data = pois, cutoff=5000)

# fit by eye!
pdf("fig/samplevarP.pdf")
plot(samplevarP)
dev.off()

#
modelvarP <- vgm(psill = 2, model = "Sph", range = 3000, nugget = 4)
modelvarP <- fit.variogram(samplevarP, model = modelvarP)
linevarP<-variogramLine(modelvarP,maxdist=max(samplevarP$dist))

# estimate hyperparameters with reml

##----bin.reml
bin.geo<-as.geodata(
    obj = as.data.frame(bin),
    header = TRUE,
    coords.col = c("x","y"),
    data.col = "residB", 
    data.names = NULL,
    covar.col = c("silt", "silt2", "depth")
)

bin.reml <- likfit(
    geodata = bin.geo,
    trend="cte", 
    cov.model="spherical", 
    ini.cov.pars=c(modelvarB[2,2], modelvarB[2,3]), 
    nugget=modelvarB[1,2],
    lik.method="REML"
)
saveRDS(bin.reml,"output/bin.reml.rds")

##----pois.reml
pois.geo<-as.geodata(
    obj = as.data.frame(pois),
    header = TRUE,
    coords.col = c("x","y"),
    data.col = "residP", 
    data.names = NULL,
    covar.col = c("silt", "silt2", "depth")
)

pois.reml <- likfit(
    geodata = pois.geo,
    trend="cte", 
    cov.model="spherical", 
    ini.cov.pars=c(modelvarP[2,2], modelvarP[2,3]), 
    nugget=modelvarP[1,2],
    lik.method="REML"
)
saveRDS(pois.reml,"output/pois.reml.rds")

# simulate S conditional on Y with Cholesky decomposition
# hyperparameters estimated with reml are set as initial

bin.mcmc.geo<-as.geodata(
    obj = as.data.frame(bin),
    header = TRUE,
    coords.col = c("x","y"),
    data.col = "macoma", 
    data.names = NULL,
    covar.col = c("silt", "silt2", "depth")
)

mcmcSetB <- mcmc.control(S.scale = 0.15, thin = 50, burn.in = 100)

glgmB <- list(
    family='binomial', 
    trend=trend.spatial(trend=trendB, geodata=as.data.frame(bin)),
    cov.model='spherical', 
    cov.pars=c(bin.reml$sigmasq, bin.reml$phi), 
    nugget =bin.reml$tausq,
    beta=bin.glm.beta)

bin.simFtilde <- glsm.mcmc(
    geodata = bin.mcmc.geo, 
    model = glgmB, 
    mcmc.input = mcmcSetB,
    messages=TRUE)

saveRDS(bin.simFtilde,"output/bin.simFtilde.rds")

#bin.simFtilde<-readRDS("output/bin.simFtilde.rds")
#dat.bin<-bin.simFtilde$geodata
#saveRDS(dat.bin,"output/dat.bin.rds")

bin.chainConv1 <- create.mcmc.coda(x = bin.simFtilde$simulations[round(runif(1, min = 1, max = nrow(bin.simFtilde$simulations)),0),], mcmc.input = mcmcSetB) 

pdf("./fig/bin.simFtilde.pdf")
par(mfrow=c(2,2))
traceplot(bin.chainConv1)
autocorr.plot(bin.chainConv1, auto.layout=FALSE)
densplot(bin.chainConv1)
geweke.plot(bin.chainConv1, auto.layout=FALSE)
dev.off()

##----pois.simFtilde
pois.mcmc.geo<-as.geodata(
    obj = as.data.frame(pois),
    header = TRUE,
    coords.col = c("x","y"),
    data.col = "macoma", 
    data.names = NULL,
    covar.col = c("silt", "silt2", "depth")
)

mcmcSetP <- mcmc.control(S.scale = 0.1, thin = 50, burn.in = 100)

glgmP <- list(
    family='poisson', 
    trend =trend.spatial(trendP, as.data.frame(pois)),
    cov.model='spherical', 
    cov.pars=c(pois.reml$sigmasq, pois.reml$phi), 
    nugget =pois.reml$tausq,
    beta=pois.glm.beta)

#constructed on kreukel
#load("input/glgmP.RData")

pois.simFtilde <- glsm.mcmc(
    geodata = pois.mcmc.geo, 
    model = glgmP, 
    mcmc.input = mcmcSetP,
    messages=TRUE)

saveRDS(pois.simFtilde,"output/pois.simFtilde.rds")

#pois.simFtilde<-readRDS("output/pois.simFtilde.rds")
#dat.pois<-pois.simFtilde$geodata
#saveRDS(dat.pois,"output/dat.pois.rds")

pois.chainConv1 <- create.mcmc.coda(x = pois.simFtilde$simulations[round(runif(1, min = 1, max = nrow(pois.simFtilde$simulations)),0),], mcmc.input = mcmcSetP) 

pdf("./fig/pois.simFtilde.pdf")
par(mfrow=c(2,2))
traceplot(pois.chainConv1)
autocorr.plot(pois.chainConv1, auto.layout=FALSE)
densplot(pois.chainConv1)
geweke.plot(pois.chainConv1, auto.layout=FALSE)
dev.off()

##----bin.mcml1
bin.mcmlPrep1 <- prepare.likfit.glsm(bin.simFtilde)  
bin.mcml1 <- likfit.glsm(
    mcmc.obj = bin.mcmlPrep1, 
    cov.model = "spherical", 
    ini.phi = bin.reml$phi, 
    nugget.rel = (bin.reml$tausq/bin.reml$sigmasq), 
    fix.nugget.rel = FALSE,
    messages = FALSE
)
saveRDS(bin.mcml1, "output/bin.mcml1.rds")

##----bin.simF
bin.simF <- glsm.mcmc(
    geodata =bin.mcmc.geo, 
    units.m = "default", 
    model = bin.mcml1,  
    mcmc.input = mcmcSetB, 
    messages=FALSE
)
saveRDS(bin.simF, "output/bin.simF.rds")

bin.chainConv2 <- create.mcmc.coda(x = bin.simF$simulations[round(runif(1, min = 1, max = nrow(bin.simF$simulations)),0),], mcmc.input = mcmcSetB)  

pdf("./fig/bin.simF.pdf")
par(mfrow=c(2,2))
traceplot(bin.chainConv2)
autocorr.plot(bin.chainConv2, auto.layout=FALSE)
densplot(bin.chainConv2)
geweke.plot(bin.chainConv2, auto.layout=FALSE)
dev.off()

##----pois.mcml1
pois.mcmlPrep1 <- prepare.likfit.glsm(pois.simFtilde)  
pois.mcml1 <- likfit.glsm(
    mcmc.obj = pois.mcmlPrep1, 
    cov.model = "spherical", 
    ini.phi = pois.reml$phi, 
    nugget.rel = (pois.reml$tausq/pois.reml$sigmasq), 
    fix.nugget.rel = FALSE,
    messages = FALSE
)
saveRDS(pois.mcml1, "output/pois.mcml1.rds")

##----pois.simF
pois.simF <- glsm.mcmc(
    geodata =pois.mcmc.geo, 
    units.m = "default", 
    model = pois.mcml1,  
    mcmc.input = mcmcSetP, 
    messages=FALSE
)
saveRDS(pois.simF, "output/pois.simF.rds")

pois.chainConv2 <- create.mcmc.coda(x = pois.simF$simulations[round(runif(1, min = 1, max = nrow(pois.simF$simulations)),0),], mcmc.input = mcmcSetP) 

pdf("./fig/pois.simF.pdf")
par(mfrow=c(2,2))
traceplot(pois.chainConv2)
autocorr.plot(pois.chainConv2, auto.layout=FALSE)
densplot(pois.chainConv2)
geweke.plot(pois.chainConv2, auto.layout=FALSE)
dev.off()

##----bin.mcml2
bin.mcmlPrep2 <- prepare.likfit.glsm(bin.simF)   
bin.mcml2 <- likfit.glsm(
  mcmc.obj = bin.mcmlPrep2, 
  cov.model = "spherical", 
  ini.phi = bin.mcml1$cov.pars[2], 
  nugget.rel = bin.mcml1$nugget.rel, 
  fix.nugget.rel = FALSE,
  messages = TRUE
)
saveRDS(bin.mcml2,"output/bin.mcml2.rds")
##----bin.simF2
bin.simF2 <- glsm.mcmc(
  geodata = bin.mcmc.geo, 
  units.m = "default", 
  model = bin.mcml2, 
  mcmc.input = mcmcSetB, #input parameters for MCMC
  messages=TRUE
)
saveRDS(bin.simF2,"output/bin.simF2.rds")

##----pois.mcml2
pois.mcmlPrep2 <- prepare.likfit.glsm(pois.simF)   
pois.mcml2 <- likfit.glsm(
  mcmc.obj = pois.mcmlPrep2, 
  cov.model = "spherical", 
  ini.phi = pois.mcml1$cov.pars[2], 
  nugget.rel = pois.mcml1$nugget.rel, 
  fix.nugget.rel = FALSE,
  messages = TRUE
)
saveRDS(pois.mcml2,"output/pois.mcml2.rds")

##----pois.simF2
pois.simF2 <- glsm.mcmc(
  geodata = pois.mcmc.geo, 
  units.m = "default", 
  model = pois.mcml2, 
  mcmc.input = mcmc.control(S.scale = 0.2, thin = 100, burn.in = 100),
  messages=TRUE
)
saveRDS(pois.simF2,"output/pois.simF2.rds")

##----pois.mcml2
#extra step for poisson only
pois.mcmlPrep3 <- prepare.likfit.glsm(pois.simF2)   
pois.mcml3 <- likfit.glsm(
  mcmc.obj = pois.mcmlPrep3, 
  cov.model = "spherical", 
  ini.phi = pois.mcml2$cov.pars[2], 
  nugget.rel = pois.mcml2$nugget.rel, 
  fix.nugget.rel = FALSE,
  messages = TRUE
)
saveRDS(pois.mcml3,"output/pois.mcml3.rds")

##----pois.simF3
pois.simF3 <- glsm.mcmc(
  geodata = pois.mcmc.geo, 
  units.m = "default", 
  model = pois.mcml3, 
  mcmc.input = mcmc.control(S.scale = 0.2, thin = 100, burn.in = 100),
  messages=TRUE
)
saveRDS(pois.simF3,"output/pois.simF3.rds")

# make 100 x 100 prediction grid 

##----sgrid
library(rgdal)
library(rgeos)
library(raster)
library(maptools)

#projection
RD<-CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")

#load shape file of study area
intertidal.shp<-readOGR("input/intertidal", "intertidal")
intertidal.shp<-spTransform(intertidal.shp, RD)
#writePolyShape(intertidal.shp, "tmp")
sgrid <- spsample(
      intertidal.shp,
      cellsize =100,
      type = "regular",
      offset = c(0.5, 0.5)) #centric systematic
proj4string(sgrid) <-RD
depth<-raster("input/depth2.grd")
mindepth<- -200
maxdepth<- 100
grddepth <- over(sgrid, as(depth,"SpatialGridDataFrame"))
ids<-which(grddepth$layer> mindepth & grddepth$layer < maxdepth )
sgrid<-sgrid[ids,]
grddepth<-grddepth[ids,]
#mgs
tmp<-dat
coordinates(tmp)<-~x+y
proj4string(tmp)<-RD
#idw mgs scaled, no need to scale again
silt.idw<-idw(silt~1, tmp, sgrid)
silt2.idw<-idw(silt2~1, tmp, sgrid)
sgrid<-as.data.frame(sgrid)
sgrid$silt<-silt.idw$var1.pred
sgrid$silt2<-silt2.idw$var1.pred
sgrid$depth<-(grddepth-meancov["depth"])/sdcov["depth"]

sgrid<-sgrid[complete.cases(sgrid),]
coordinates(sgrid)<-c("x","y")
proj4string(sgrid) <- RD
saveRDS(sgrid, "output/sgrid.rds")

