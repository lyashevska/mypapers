# select one pseudo-reality that resembles raw data SIBES using the following criteria:
# load simulated signals on 100 m grid and 1000 validation points
# percentage zeros
# mean of Poisson part
# var of poisson part

rm(list=ls())
set.seed(5021)

# load packages, functions
library(geoRglm)
library(gstat)
library(sp)
library(ggplot2)
library(gstat)
library(pscl)
library(rgdal)
library(reshape)

# projection rijksdriehook
RD<-CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")

# load custom functions
source("./myfun/myfun.count.R")
#source("./myfun/myfun.sample.R")

# use antilogit to transform simulated S to Binomial pi and exp to transform signals to Poisson mu
antilogit<-function(x){exp(x)/(1+exp(x))} #for bin

# use to transford kriging interpolation (prediction) to pi and mu
transf.predP<-function(var1.pred, var1.var){exp(var1.pred+0.5*var1.var)}
transf.predB<-function(var1.pred, var1.var){plogis(var1.pred) + 0.5*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)*var1.var}

transf.varP<-function(var1.pred, var1.var){(exp(2*var1.pred+var1.var))*expm1(var1.var)}
transf.varB<-function(var1.pred, var1.var){(exp(var1.pred)/(1+exp(var1.pred))^2)^2*var1.var+(1/2)*(exp(var1.pred)*(-expm1(var1.pred))/(1+exp(var1.pred))^3)^2*var1.var^2}

## select X on 100 m that best resembles sibes
## function to define criteria to select simulated field that resembles sibes
##criteria 1: fraction of zeros
#f1<-function(x){table(x>1)[1]/length(x)}
##criteria 2: mean of non zero
#f2<-function(x){mean(subset(x,x!=0))}
##criteria 3: var of non zero
#f3<-function(x){var(subset(x,x!=0))}
#
## load scaled sibes data
#dat<-read.csv(file="./input/datsc.csv")
#r1<-apply(dat[2],2,f1)
##0.80
#r2<-apply(dat[2],2,f2)
##4.13
#r3<-apply(dat[2],2,f3)
##59.24
#
## load signals on 100 m grid 
#simB<-readRDS("./output/simB.rds")
#simP<-readRDS("./output/simP.rds")
#
## transform simulation
## exclude coordinate columns, first two
#sim.pi<-antilogit(simB[,-c(1:2)])
#sim.mu<-exp(simP[,-c(1:2)])
## keep coordinate columns
#xy<-simB[,c(1:2)]
#
## obtain 100 field with counts
#sim.count<-myfun.count(sim.pi,sim.mu)
## apply criteria
#r1X<-apply(sim.count,2,f1)
#r2X<-apply(sim.count,2,f2)
#r3X<-apply(sim.count,2,f3)
## divide abs difference by standard deviation
#r1Xa<- abs(r2-r1X)/(sd(abs(r1-r1X)))
#r2Xa<- abs(r2-r2X)/(sd(abs(r2-r2X)))
#r3Xa<- abs(r3-r3X)/(sd(abs(r3-r3X)))
## find sum of absolute differences
## assign weights 0.5, 0.25, 0.25
## and select min
#rX<-0.5*r1Xa+0.25*r2Xa+0.25*r3Xa
#which.min(rX)[[1]]
##summary(rX)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  541.3   543.0   543.6   543.6   544.1   545.9 
## simulated pseudoreality X69 mostly resembles SIBES 
## with fraction of zeros being 0.8158977, 
## mean of non-zero 4.131529 
## variance of non-zero 37.47579 
#
##realisation 69 resembles raw SIBES 
##select X69
#
###prediction grid with covariates
##sgrid<-readRDS("input/sgrid.rds")
###set aside 1000 points based on tidal basins
###load tidal basins
##tidalbasins.shp<-readOGR("./input/tidalbasins", "tidalbasins")
###reproject tidal basin
##tidalbasins.shp<-spTransform(tidalbasins.shp, RD)
##tb<-over(sgrid,tidalbasins.shp)
##sgriddf<-as.data.frame(sgrid)
##sgriddf$tb<-tb$TB
###compute size (number of points) of strata
##Ntb<-table(tb$TB)
##N<-sum(Ntb)
###compute sample sizes (proportional allocation)
##n<-1000
##ntb<-round(Ntb/N*n)
##sum(ntb)
###in stratum 1 only 1 point, change this to 2.
##ids<-which(ntb<2)
##ntb[ids]<-2
##sumn<-sum(ntb)
##ids<-which(ntb==max(ntb))
##ntb[ids]<-ntb[ids]-(sum(ntb)-n)
##sum(ntb)
##stratumid<-sort(unique(tb$TB))
##val<-NULL
###set aside 1000 points using tidal basins as strata
##for (i in stratumid){
##sam<-sgriddf[myfun.sample(which(sgriddf$tb==i),ntb[[as.character(i)]], replace=FALSE),] 
##val<-rbind(val, sam)
##}
###val.id
##sgriddf$val.id<-row.names(sgriddf)%in%row.names(val)
#
##saveRDS(sgriddf, "output/sgriddf.rds")
#
#sgriddf<-readRDS("output/sgriddf.rds")
#val.xy<-subset(sgriddf, val.id==T, select=c(x1,x2))
#
##add covariates
#simBdf<-data.frame(S=simB$S69,sgriddf)
#simPdf<-data.frame(S=simP$S69,sgriddf)
#
##prediction dataset
#grdB<-simBdf[which(simBdf$val.id==F),]
#grdP<-simPdf[which(simPdf$val.id==F),]
#
##validation dataset
#valB<-simBdf[which(simBdf$val.id==T),]
#valP<-simPdf[which(simPdf$val.id==T),]
#
##assign coordinates
#coordinates(valB)<-~x1+x2
#proj4string(valB)<-RD
#coordinates(valP)<-~x1+x2
#proj4string(valP)<-RD
#coordinates(grdB)<-~x1+x2
#proj4string(grdB)<-RD
#coordinates(grdP)<-~x1+x2
#proj4string(grdP)<-RD
#gridded(grdB)<-TRUE
#gridded(grdP)<-TRUE
#
#saveRDS(grdB, file="output/grdB.rds")
#saveRDS(grdP, file="output/grdP.rds")
#saveRDS(valB, file="output/valB.rds")
#saveRDS(valP, file="output/valP.rds")
#
################ RUN HERE

# load simulated signals on 100 m grid and 1000 validation points
grdB<- readRDS("output/grdB.rds")
grdP<- readRDS("output/grdP.rds")
valB<- readRDS("output/valB.rds")
valP<- readRDS("output/valP.rds")

#MCML parameters evaluated on 100m grid
mcmlEstimation2B<-readRDS("output/bin.mcml2.rds")
mcmlEstimation2P<-readRDS("output/pois.mcml3.rds")

# define trend and model
trendB<-trendP<- S ~ silt + silt2 + depth
vgmB <- vgm(model="Sph",
  psill = mcmlEstimation2B$cov.pars[1],
  range= mcmlEstimation2B$cov.pars[2],
  nugget= mcmlEstimation2B$nugget.rel*mcmlEstimation2B$cov.pars[1])
vgmP <- vgm(model="Sph",
  psill = mcmlEstimation2P$cov.pars[1],
  range= mcmlEstimation2P$cov.pars[2],
  nugget= mcmlEstimation2P$nugget.rel*mcmlEstimation2P$cov.pars[1])

#grid spacing
spacing<-2^seq(1,5,1)*100 
# steekrproef number
nsteek<-100

seB<-seP<-se.pi<-se.mu<-prB<-prP<-pr.pi<-pr.mu<-var.pi<-var.mu<-array(dim=c(length(spacing), nsteek, 1000))
npointsB<-npointsP<-array(dim=c(length(spacing),nsteek))

#binomial
for(k in 1:length(spacing)){
for (j in 1:nsteek){
s<- spsample(grdB, cellsize=spacing[k],type="regular")
npt<-length(s)
o<-over(s, grdB)
sam<-data.frame(s,o)
names(sam)[1:2]<-c("x1","x2")
sam<-na.omit(sam)
coordinates(sam)<-~x1+x2
proj4string(sam)<-RD
dataB<-sam
tmpB <- krige(
  formula = trendB,
  locations = dataB,
  newdata=valB,
  model = vgmB,
  beta = mcmlEstimation2B$beta,
  nmax=100,
  debug.level=1
)
#signal
prB[k,j,]<-tmpB@data$var1.pred
seB[k,j,]<-(tmpB@data$var1.pred-valB$S)^2
#pi
pr.pi[k,j,]<-transf.predB(tmpB@data$var1.pred,tmpB@data$var1.var)
#kriging variance
var.pi[k,j,]<-transf.varB(tmpB@data$var1.pred,tmpB@data$var1.var)
se.pi[k,j,]<-(transf.predB(tmpB@data$var1.pred, tmpB@data$var1.var)-antilogit(valB$S))^2
npointsB[k,j]<-npt
}
}

#poisson
for(k in 1:length(spacing)){
for (j in 1:nsteek){
s<- spsample(grdP, cellsize=spacing[k],type="regular")
npt<-length(s)
o<-over(s, grdP)
sam<-data.frame(s,o)
names(sam)[1:2]<-c("x1","x2")
sam<-na.omit(sam)
coordinates(sam)<-~x1+x2
proj4string(sam)<-RD
dataP<-sam
tmpP <- krige(
  formula = trendP,
  locations = dataP,
  newdata=valP,
  model = vgmP,
  beta = mcmlEstimation2P$beta,
  nmax=100,
  debug.level=1
)
#signal
prP[k,j,]<-tmpP@data$var1.pred
seP[k,j,]<-(tmpP@data$var1.pred-valP$S)^2
#mu
pr.mu[k,j,]<-transf.predP(tmpP@data$var1.pred,tmpP@data$var1.var)
#kriging variance
var.mu[k,j,]<-transf.varP(tmpP@data$var1.pred,tmpP@data$var1.var)
se.mu[k,j,]<-(transf.predP(tmpP@data$var1.pred, tmpP@data$var1.var)-exp(valP$S))^2
npointsP[k,j]<-npt
}
}

#calculate mse for signals B/P and for pi/mu
mseB<-apply(seB,c(1:2),mean)
mse.pi<-apply(se.pi,c(1:2),mean)
mseP<-apply(seP,c(1:2),mean)
mse.mu<-apply(se.mu,c(1:2),mean) # see below for bernoulli 1 only!
#1000 val points averaged over nsteek per grid spacing
mean.pr.pi<-apply(pr.pi,c(1,3),mean)
mean.pr.mu<-apply(pr.mu,c(1,3),mean)
val.pi<- antilogit(valB$S)
save(seB,seP,prB,prP, pr.pi, pr.mu,se.pi,se.mu,var.pi,var.mu, mse.pi, mse.mu, npointsB, npointsP,val.pi,
     file="./output/results.RData")

##calculate mse.mu only for validation points with Bernoulli 1
val.pi<-antilogit(valB$S)
#vector where bernoulli 1 = TRUE 
val.ber1<-val.pi>0.5
saveRDS(val.ber1, "output/val.ber1.rds")

nsteek<-100
se.mu1<-array(dim=c(length(spacing), nsteek, length(val.ber1[val.ber1==TRUE])))

for(k in 1:length(spacing)){
for (j in 1:nsteek){
se.mu1[k,j,]<-se.mu[k,j,][val.ber1]
}
}
mse.mu1<-apply(se.mu1,c(1:2),mean)
saveRDS(mse.mu1, "output/mse.mu1.rds")
