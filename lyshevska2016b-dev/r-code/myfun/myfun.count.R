myfun.count<-function(sim.pi,sim.mu){
datb<-numeric()
for (i in 1:ncol(sim.pi)){
#generate random samples from bin dist
b<-rbinom(n=nrow(sim.pi), size=1, p=sim.pi[,i])
datb<-cbind(datb, b)}
datp<-numeric()
for (j in 1:ncol(sim.mu)){
#generate random samples from pois dist
d<-rpois(n=nrow(sim.mu), lambda=sim.mu[,j])
datp<-cbind(datp, d)}
res<-ifelse(datb==0,0,datp)
res<-data.frame(res)
names(res)<- c(paste(rep("X",ncol(sim.pi)),c(1:ncol(sim.pi)),sep=""))
res<-res[complete.cases(res),]
}

