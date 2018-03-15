#function to partition data
#separate zeros: bernoulli and poisson
myfun.sep<-function(sam){
samB<-sam
samB$macoma<-ifelse(samB$macoma>0,1,0)
samP<-sam
glmB <- glm(
    formula = macoma~silt+silt2+depth,
    family = "binomial", 
    data = subset(samB, select=-c(x,y))
)
glmZ <- zeroinfl(
    formula =macoma~silt+silt2+depth,
    link = "logit",
    dist="poisson",
    data = samP
)
p.glmB <- predict(glmB, type='response') 
summary(p.glmB)
lp.glmZ <- predict(glmZ, type='zero') 
summary(lp.glmZ)
p.glmZ<-1-lp.glmZ 
summary(p.glmZ)
mu.glmZ <- predict(glmZ,type='count') 
p.glmZ.tot <- p.glmZ*exp(-mu.glmZ)
p.rat <- (1-p.glmZ)/(1-p.glmZ+p.glmZ.tot)
p.test <- runif(length(samP$macoma))
echte <- p.rat > p.test
keep<-ifelse(echte==TRUE & sam$macoma==0, FALSE, TRUE)
samP0<-samP[keep,]
samB0<-samB
samB0$macoma<-ifelse(keep==TRUE, 1, samB$macoma)
return(list(bin=samB0,pois=samP0))
}
