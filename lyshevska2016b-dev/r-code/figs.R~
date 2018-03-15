###### FIGURES
# AbundanceRaw
# pred.pimu
# mse.mu1
# mse.pe
# mse.pimu

# load libraries
rm(list=ls())
library(ggplot2)
library(reshape2)
library(grid)
library(sp)
library(rgdal)

# grid spacing
spacing<-2^seq(1,5,1)*100 

# projection rijksdriehook
RD<-CRS("+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812 +units=m +no_defs")
antilogit<-function(x){exp(x)/(1+exp(x))} #for bin


#load R1 MCML results
load(file="./output/R1/results.RData")
mse.mu.mcml<-mse.mu[1:5,]
mse.pi.mcml<-mse.pi[1:5,]

# average over 100
mse.pi.mcml.avg<-rowMeans(mse.pi.mcml)
mse.pimu.mcml.avg<-rowMeans(mse.pi.mcml*mse.mu.mcml)

load(file="output/R2/results.RData")
mse.mu.reml1<-apply(se.mu,1,mean)
mse.pi.reml1<-apply(se.pi,1,mean)
#rm(list= ls()[!(ls() %in% c('mse.mu.reml1','mse.pi.reml1'))])
pr.pi.reml1<-pr.pi
pr.mu.reml1<-pr.mu

load(file="output/R3/results.RData")
mse.mu.reml2<-apply(se.mu,1,mean)
mse.pi.reml2<-apply(se.pi,1,mean)
pr.pi.reml2<-pr.pi
pr.mu.reml2<-pr.mu

load(file="output/R4/results.RData")
mse.mu.reml3<-apply(se.mu,1,mean)
mse.pi.reml3<-apply(se.pi,1,mean)
pr.pi.reml3<-pr.pi
pr.mu.reml3<-pr.mu

load(file="output/R5/results.RData")
mse.mu.reml4<-apply(se.mu,1,mean)
mse.pi.reml4<-apply(se.pi,1,mean)
pr.pi.reml4<-pr.pi
pr.mu.reml4<-pr.mu

load(file="output/R6/results.RData")
mse.mu.reml5<-apply(se.mu,1,mean)
mse.pi.reml5<-apply(se.pi,1,mean)
pr.pi.reml5<-pr.pi
pr.mu.reml5<-pr.mu

# load(file="output/R6/results.RData")

#rm(list= ls()[!(ls() %in% c('mse.mu.reml1','mse.pi.reml1', 'mse.mu.reml2','mse.pi.reml2', 'mse.mu.reml3','mse.pi.reml3','mse.mu.reml4','mse.pi.reml4'))])

# calculate avg over 5 reml runs, for papers
mse.pi.reml.avg<-colMeans(rbind(mse.pi.reml1,mse.pi.reml2, mse.pi.reml3, mse.pi.reml4, mse.pi.reml5))
mse.pimu.reml.avg<-colMeans(rbind(mse.pi.reml1*mse.mu.reml1,mse.pi.reml2*mse.mu.reml2, mse.pi.reml3*mse.mu.reml3, mse.pi.reml4*mse.mu.reml4, mse.pi.reml5*mse.mu.reml5))

### gg.intertidal
intertidal.shp<-readOGR("./input/intertidal", "intertidal")
intertidal.shp<-spTransform(intertidal.shp, RD)
gg.intertidal<-ggplot(intertidal.shp) + 
geom_polygon(aes(x = long/1000, y = lat/1000, group=group), alpha=0.3)+
scale_y_continuous(name = "Northing (km)") + scale_x_continuous(name = "Easting (km)") + coord_equal()

###figure AbundanceRaw (1a)
df<-read.csv("input/dat.csv")
b<-c(-Inf,1,2, Inf)
l<-c("[0-1)","[1-2)", "[2-84)")
# colo<-c( "#FFFF00", "#7FC400", "#008B00")
# bpy.colors(n=3)
#colo<-c("#FFFF60FF", "#C729D6FF", "#000033FF")
colo<-c("#FFFFCC" , "#CC6600", "#000033FF")

#breaks for 
df$group<-cut(df$macoma,breaks=b, labels=l)

png("fig/AbundanceRaw.png", width=12, height=6, units = "in", res=300)
gg.intertidal+
geom_point(data = df, aes( x = x* 1.0e-3, y = y* 1.0e-3, colour=group),size=0.5)+ 
scale_colour_manual(name = "", values = colo)+
theme(legend.position="bottom", panel.background = element_rect(fill="white", color=NA), 
     panel.border = element_rect(colour = "grey50", fill=NA))+ 
guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

###figure pred.pimu (1b)
#realisation of selected pseudo-reality (69) at 100m 
#          0%          25%          50%          75%         100% 
#9.839988e-04 1.482087e-01 4.088860e-01 1.075328e+00 1.611503e+02 

b<-c(-Inf,0.2,0.5,1.0, Inf)
l<-c("[0-0.2)", "[0.2-0.5)","[0.5-1.0)", "[1.0-Inf)")

# colo<-c("#FFFF00", "#AAD800", "#55B100", "#008B00")
#colo<-c("#FFFF60FF", "#FF758AFF", "#5000FFFF", "#000033FF")

colo<-c("#FFFFCC" , "#FFFF60FF", "#CC6600", "#000033FF")

#sgriddf<-readRDS("output/sgriddf.rds")

# load signals on 100 m grid 
simB<-readRDS("./output/simB.rds")
simP<-readRDS("./output/simP.rds")

df<-data.frame(simB[, c(1:2)], pimu=antilogit(simB[,69])*exp(simP[,69]))

#use the same splits as for 100 m for comparison
df$group<-cut(df$pimu, breaks=b, labels=l, right=F) 

 
png("fig/pred.pimu.png", width=12, height=6, units = "in", res=300) 
gg.intertidal+ geom_point(data = df, aes( x = x1* 1.0e-3, y = x2* 1.0e-3, colour=group),size=0.001)+
scale_colour_manual(name="", values = colo)+
coord_equal(ratio = 1)+
theme(legend.position="bottom", panel.background = element_rect(fill="white", color=NA), 
     panel.border = element_rect(colour = "grey50", fill=NA))+ 
guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()


###figure mse.pi (2a)
row.names(mse.pi.mcml)<-spacing
df1<-melt(mse.pi.mcml)
df2<-data.frame(spacing, mse.pi.reml=c(NA, mse.pi.reml1))
df3<-data.frame(spacing, mse.pi.reml=c(NA, mse.pi.reml2))
df4<-data.frame(spacing, mse.pi.reml=c(NA, mse.pi.reml3))
df5<-data.frame(spacing, mse.pi.reml=c(NA, mse.pi.reml4))
df6<-data.frame(spacing, mse.pi.reml=c(NA, mse.pi.reml5))
png("fig/mse.pi.png", width=700*3.125, height=700*3.125, res=300) 
ggplot(df1, aes(x=as.factor(Var1), y=value))+
geom_boxplot(outlier.size=2, notch=F, fill="skyblue")+
geom_point(data=df2, aes(x=as.factor(spacing), y=mse.pi.reml), size=3, shape=8, col="red")+
geom_point(data=df3, aes(x=as.factor(spacing), y=mse.pi.reml), size=3, shape=8, col="red")+
geom_point(data=df4, aes(x=as.factor(spacing), y=mse.pi.reml), size=3, shape=8, col="red")+
geom_point(data=df5, aes(x=as.factor(spacing), y=mse.pi.reml), size=3, shape=8, col="red")+
geom_point(data=df6, aes(x=as.factor(spacing), y=mse.pi.reml), size=3, shape=8, col="red")+
scale_x_discrete(name = "Spacing (m)") + 
scale_y_continuous(name = "MSE (prevalence)")+theme_bw()
#scale_y_continuous(limits=c(0,0.02), name = "MSE (prevalence)")+theme_bw()
dev.off()

###figure mse.pimu (2c)
row.names(mse.pi.mcml)<-spacing
df1<-melt(mse.pi.mcml*mse.mu.mcml)
df2<-data.frame(spacing, mse.pimu.reml=c(NA, mse.pi.reml1*mse.mu.reml1)) 
df3<-data.frame(spacing, mse.pimu.reml=c(NA, mse.pi.reml2*mse.mu.reml2)) 
df4<-data.frame(spacing, mse.pimu.reml=c(NA, mse.pi.reml3*mse.mu.reml3)) 
df5<-data.frame(spacing, mse.pimu.reml=c(NA, mse.pi.reml4*mse.mu.reml4)) 
df6<-data.frame(spacing, mse.pimu.reml=c(NA, mse.pi.reml5*mse.mu.reml5)) 
png("fig/mse.pimu.png", width=700*3.125, height=700*3.125, res=300) 
ggplot(df1, aes(x=as.factor(Var1), y=value))+geom_boxplot(outlier.size=2, notch=F, fill="skyblue")+
geom_point(data=df2, aes(x=as.factor(spacing), y=mse.pimu.reml), size=3, shape=8, col="red")+
geom_point(data=df3, aes(x=as.factor(spacing), y=mse.pimu.reml), size=3, shape=8, col="red")+
geom_point(data=df4, aes(x=as.factor(spacing), y=mse.pimu.reml), size=3, shape=8, col="red")+
geom_point(data=df5, aes(x=as.factor(spacing), y=mse.pimu.reml), size=3, shape=8, col="red")+
geom_point(data=df6, aes(x=as.factor(spacing), y=mse.pimu.reml), size=3, shape=8, col="red")+
scale_x_discrete(name = "Spacing (m)") + 
scale_y_continuous(name = "MSE (unconditional intensity)")+
#scale_y_continuous(limits=c(0,0.1), name = "MSE (unconditional intensity)")+
theme_bw()
dev.off()

###figure mse.mu1 (2b)
#mse.mu1 graph, with Bernoulli 1
mse.mu<-readRDS("./output/R1/mse.mu1.rds")
mse.mu.mcml<-mse.mu[1:5,]
#load validation set for bernoulli 1
val.ber1<-readRDS("./output/R1/val.ber1.rds")
#load reml results
load(file="output/R2/results.RData")
se.mu1<-array(dim=c(length(spacing), length(val.ber1[val.ber1==TRUE])))
for(k in 1:(length(spacing)-1)){
se.mu1[k,]<-se.mu[k,][val.ber1]
}
mse.mu.reml1<-apply(se.mu1,1,mean)[-c(5,6)]

##load reml results
load(file="output/R3/results.RData")
se.mu1<-array(dim=c(length(spacing), length(val.ber1[val.ber1==TRUE])))
for(k in  1:(length(spacing)-1)){
se.mu1[k,]<-se.mu[k,][val.ber1]
}
mse.mu.reml2<-apply(se.mu1,1,mean)[-c(5,6)]

#load reml results
load(file="output/R4/results.RData")
se.mu1<-array(dim=c(length(spacing), length(val.ber1[val.ber1==TRUE])))
for(k in  1:(length(spacing)-1)){
se.mu1[k,]<-se.mu[k,][val.ber1]
}
mse.mu.reml3<-apply(se.mu1,1,mean)[-c(5,6)]

#load reml results
load(file="output/R5/results.RData")
se.mu1<-array(dim=c(length(spacing), length(val.ber1[val.ber1==TRUE])))
for(k in  1:(length(spacing)-1)){
se.mu1[k,]<-se.mu[k,][val.ber1]
}
mse.mu.reml4<-apply(se.mu1,1,mean)[-c(5,6)]

#load reml results
load(file="output/R6/results.RData")
se.mu1<-array(dim=c(length(spacing), length(val.ber1[val.ber1==TRUE])))
for(k in  1:(length(spacing)-1)){
se.mu1[k,]<-se.mu[k,][val.ber1]
}
mse.mu.reml5<-apply(se.mu1,1,mean)[-c(5,6)]

#measure of quality for mu (mse mu) 
row.names(mse.mu.mcml)<-spacing
df1<-melt(mse.mu.mcml)
df2<-data.frame(spacing, mse.mu.reml=c(NA, mse.mu.reml1))
df3<-data.frame(spacing, mse.mu.reml=c(NA, mse.mu.reml2))
df4<-data.frame(spacing, mse.mu.reml=c(NA, mse.mu.reml3))
df5<-data.frame(spacing, mse.mu.reml=c(NA, mse.mu.reml4))
df6<-data.frame(spacing, mse.mu.reml=c(NA, mse.mu.reml5))

png("fig/mse.mu1.png", width=700*3.125, height=700*3.125, res=300) 
ggplot(df1, aes(x=as.factor(Var1), y=value))+
geom_boxplot(outlier.size=2, notch=F, fill="skyblue")+
geom_point(data=df2, aes(x=as.factor(spacing), y=mse.mu.reml), size=3, shape=8, col="red")+
geom_point(data=df3, aes(x=as.factor(spacing), y=mse.mu.reml), size=3, shape=8, col="red")+
geom_point(data=df4, aes(x=as.factor(spacing), y=mse.mu.reml), size=3, shape=8, col="red")+
geom_point(data=df5, aes(x=as.factor(spacing), y=mse.mu.reml), size=3, shape=8, col="red")+
geom_point(data=df6, aes(x=as.factor(spacing), y=mse.mu.reml), size=3, shape=8, col="red")+
scale_x_discrete(name = "Spacing (m)") + 
scale_y_continuous(name = "MSE (intensity)")+theme_bw()
#scale_y_continuous(limits=c(0,22),name = "MSE (intensity)")+theme_bw()
dev.off()

mse.mu.mcml.avg <- rowMeans(mse.mu.mcml)
mse.mu.reml.avg<-colMeans(rbind(mse.mu.reml1,mse.mu.reml2, mse.mu.reml3, mse.mu.reml4))

# save numbers for the papers
save(mse.pi.mcml.avg, mse.pi.reml.avg, mse.pimu.mcml.avg, mse.pimu.reml.avg, mse.mu.mcml.avg, mse.mu.reml.avg, file="output/all.mse.avg.RData")
