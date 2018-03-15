#reml
#specific parameters at each grid spacing

set.seed(22801)

#library
library(sp)
library(sp)
library(geoRglm)
library(gstat)
library(rgdal)
library(reshape)

#transformation from S to pi and mu according to christensen
transf.predB<-function(var1.pred, var1.var){plogis(var1.pred) + 0.5*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)*var1.var}
transf.varB<-function(var1.pred, var1.var){
        (exp(var1.pred)/(1+exp(var1.pred))^2)^2*var1.var+(1/2)*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)^2*var1.var^2}
transf.predP<-function(var1.pred, var1.var){exp(var1.pred+0.5*var1.var)}
transf.varP<-function(var1.pred, var1.var){(exp(2*var1.pred+var1.var))*expm1(var1.var)}

antilogit<-function(x){exp(x)/(1+exp(x))} #for bin


#input files
RD<-CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")

# load simulated signals on 100 m grid and 1000 validation points
grdB<-readRDS("output/grdB.rds")
valB<-readRDS("output/valB.rds")
grdP<-readRDS("output/grdP.rds")
valP<-readRDS("output/valP.rds")

#coords for validation points
val.xy<-data.frame(valB@coords)
#back-transform signals B and P into parameters pi and mu
grd.pi<-grdB
grd.pi$S<-antilogit(grdB$S)
grd.mu<-grdP
grd.mu$S<-exp(grdP$S)
val.pi<-valB
val.pi$S<-antilogit(valB$S)
val.mu<-valP
val.mu$S<-antilogit(valP$S)

#MCML parameters evaluated on 100m grid
mcmlEstimation2B<-readRDS("output/bin.mcml2.rds")
mcmlEstimation2P<-readRDS("output/pois.mcml3.rds")

#grid spacing
#start with 400m
#HERE do only large grid<--HERE
spacing<-2^seq(2,5,1)*100 
#spacing<-2^seq(4,5,1)*100 
k<-1

#define trend for kriging, likfit is the same but without S
trendB<-trendP<- S ~ silt + silt2 + depth

seB<-seP<-se.pi<-se.mu<-prB<-prP<-pr.pi<-pr.mu<-var.pi<-var.mu<-array(dim=c(length(spacing), 1000))
remlB.conv<-remlP.conv<-remlB.sigmasq<-remlP.sigmasq<-remlB.nugget<-remlP.nugget<-remlB.phi<-remlP.phi<-npointsB<-npointsP<-array(dim=c(length(spacing)))
remlB.beta<-remlP.beta<-array(dim=c(length(spacing),4)) 

for(k in 1:length(spacing)){
#binomial
s<- spsample(grdB, cellsize=spacing[k],type="regular")
nptB<-length(s)
o<-over(s, grdB)
sam<-data.frame(s,o)
names(sam)[1:2]<-c("x1","x2")
sam<-na.omit(sam)
coordinates(sam)<-~x1+x2
proj4string(sam)<-RD
dataB<-sam
#poisson
s<- spsample(grdP, cellsize=spacing[k],type="regular")
nptP<-length(s)
o<-over(s, grdP)
sam<-data.frame(s,o)
names(sam)[1:2]<-c("x1","x2")
sam<-na.omit(sam)
coordinates(sam)<-~x1+x2
proj4string(sam)<-RD
dataP<-sam
#estimate parameters with reml
remlB <- likfit(
    geodata = as.geodata(sam,data.col="S",covar.col=c("silt", "silt2", "depth")),
    trend=trend.spatial(~ silt + silt2 + depth), 
    cov.model="spherical", 
    ini.cov.pars=c(mcmlEstimation2B$cov.pars[1], mcmlEstimation2B$cov.pars[2]), 
    nugget=mcmlEstimation2B$nugget.rel*mcmlEstimation2B$cov.pars[1],
    lik.method="REML"
)
remlP <- likfit(
    geodata = as.geodata(sam,data.col="S",covar.col=c("silt", "silt2", "depth")),
    trend=trend.spatial(~ silt + silt2 + depth), 
    cov.model="spherical", 
    ini.cov.pars=c(mcmlEstimation2P$cov.pars[1], mcmlEstimation2P$cov.pars[2]), 
    nugget=mcmlEstimation2P$nugget.rel*mcmlEstimation2P$cov.pars[1],
    lik.method="REML"
)
#check convergence
if(remlB$info.minimisation.function$convergence==0&remlP$info.minimisation.function$convergence==0){
if(remlB$phi==0) remlB$phi<-0.001
if(remlP$phi==0) remlP$phi<-0.001
#kriging S Bernoulli
tmpB <- krige(
  formula = trendB,
  locations = dataB,
  newdata=valB,
  model = vgm(model="Sph",psill=remlB$sigmasq,range=remlB$phi,nugget=remlB$nugget),
  beta = remlB$beta,
  nmax=100,
  debug.level=1
)
#signal
prB[k,]<-tmpB@data$var1.pred
seB[k,]<-(tmpB@data$var1.pred-valB$S)^2
#pi
pr.pi[k,]<-transf.predB(tmpB@data$var1.pred,tmpB@data$var1.var)
#kriging variance
var.pi[k,]<-transf.varB(tmpB@data$var1.pred,tmpB@data$var1.var)
se.pi[k,]<-(transf.predB(tmpB@data$var1.pred, tmpB@data$var1.var)-antilogit(valB$S))^2
npointsB[k]<-nptB
#kriging S Poisson
tmpP <- krige(
  formula = trendP,
  locations = dataP,
  newdata=valP,
  model = vgm(model="Sph",psill=remlP$sigmasq,range=remlP$phi,nugget=remlP$nugget),
  beta = remlP$beta,
  nmax=100,
  debug.level=1
)
#signal
prP[k,]<-tmpP@data$var1.pred
seP[k,]<-(tmpP@data$var1.pred-valP$S)^2
#pi
pr.mu[k,]<-transf.predP(tmpP@data$var1.pred,tmpP@data$var1.var)
#kriging variance
var.mu[k,]<-transf.varP(tmpP@data$var1.pred,tmpP@data$var1.var)
se.mu[k,]<-(transf.predP(tmpP@data$var1.pred, tmpP@data$var1.var)-exp(valP$S))^2
npointsP[k]<-nptP
#accumulate reml convergence
remlB.conv[k]<-remlB$info.minimisation.function$convergence
remlB.phi[k]<-remlB$phi
remlB.sigmasq[k]<-remlB$sigmasq
remlB.nugget[k]<-remlB$nugget
remlP.conv[k]<-remlP$info.minimisation.function$convergence
remlP.phi[k]<-remlP$phi
remlP.sigmasq[k]<-remlP$sigmasq
remlP.nugget[k]<-remlP$nugget
remlB.beta[k,]<-remlB$beta
remlP.beta[k,]<-remlP$beta
}
}
mse.pi<-apply(se.pi,1,mean)
mse.mu<-apply(se.mu,1,mean)
mse.pi
mse.mu

val.pi<-antilogit(valB$S)
val.mu<-exp(valP$S)

#add results400
#save all results 
save(mse.pi, mse.mu, seB,seP,prB,prP, pr.pi, pr.mu,se.pi,se.mu,var.pi,var.mu, npointsB, npointsP,val.xy, val.pi, val.mu,
     file="./output/R2/results.RData")
#save reml parameters
save(remlB.conv, remlP.conv, remlB.sigmasq, remlP.sigmasq, remlB.nugget, remlP.nugget, npointsB, npointsP, remlB.beta, remlP.beta,
     file="./output/R2/reml.RData")

