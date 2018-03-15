# conditional sequential gausssian simulation 101 m grid

rm(list=ls())
set.seed(5021)

# load libraries
library(gstat)
library(geoR)
library(raster)
library(sp)
library(geoRglm)
library(ggplot2)
library(gstat)
library(pscl)
library(MASS)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)

# load 
simF2B<- readRDS("./output/bin.simF2.rds")
simF2P<- readRDS("./output/pois.simF3.rds")
mcmlEstimation2B<-readRDS("./output/bin.mcml2.rds")
mcmlEstimation2P<-readRDS("./output/pois.mcml3.rds")
sgrid<-readRDS("./output/sgrid.rds")

# projection rijksdriehook
RD<-CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")

##----vgmB
#setup prediction model
vgmB <- vgm(model="Sph",
psill = mcmlEstimation2B$cov.pars[1],
range= mcmlEstimation2B$cov.pars[2],
nugget= mcmlEstimation2B$nugget.rel*mcmlEstimation2B$cov.pars[1])

##----vgmP
vgmP <- vgm(model="Sph",
psill = mcmlEstimation2P$cov.pars[1],
range= mcmlEstimation2P$cov.pars[2],
nugget= mcmlEstimation2P$nugget.rel*mcmlEstimation2P$cov.pars[1])

##----simB
nsm<-100
simB<-NULL
for (i in 1:nsm){
S.sibesB<-as.data.frame(cbind(simF2B$geodata$coords,simF2B$simulations[,i], simF2B$geodata$covariate))
names(S.sibesB)<-c("x","y","S", names(simF2B$geodata$covariate))
coordinates(S.sibesB)<-~x+y
proj4string(S.sibesB)<-RD
S.krige.simB <- krige(
    formula = S~silt+silt2+depth,
    locations = S.sibesB,
    newdata=sgrid,
    model = vgmB,
    beta = mcmlEstimation2B$beta,
    nmax=100,
    nsim=1, 
    debug.level=1
)

#simulation
simB<-cbind(simB,S.krige.simB@data$sim1)
}
colnames(simB)<-paste(rep("S",ncol(simB)),c(1:ncol(simB)),sep="")
simB<-data.frame(sgrid@coords, simB)
saveRDS(object=simB,file="./output/simB.rds")

##----simP
nsm<-100
simP<-NULL
for (i in 1:nsm){
S.sibesP<-as.data.frame(cbind(simF2P$geodata$coords,simF2P$simulations[,i], simF2P$geodata$covariate))
names(S.sibesP)<-c("x","y","S", names(simF2P$geodata$covariate))
coordinates(S.sibesP)<-~x+y
proj4string(S.sibesP)<-RD
S.krige.simP <- krige(
    formula = S~silt+silt2+depth,
    locations = S.sibesP,
    newdata=sgrid,
    model = vgmP,
    beta = mcmlEstimation2P$beta,
    nmax=100,
    nsim=1, 
    debug.level=1
)

#simulation
simP<-cbind(simP,S.krige.simP@data$sim1)
}
colnames(simP)<-paste(rep("S",ncol(simP)),c(1:ncol(simP)),sep="")
simP<-data.frame(sgrid@coords, simP)
saveRDS(object=simP,file="./output/simP.rds")

