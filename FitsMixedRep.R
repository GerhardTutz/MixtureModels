


### mixture model   fits 

### 
#getwd() # check 
#setwd("C:\\Users\\tutz\\LRZ Sync+Share\\ABookOrdinal\\R\\Mixtures\\")


library("ordinal")
library("VGAM")
library("xtable")
#library("ordinalgmifs")
library("ordDisp")
library("CUB")
install.packages("FastCUB")
library("FastCUB")
library("ordDisp")
library(MASS)


### use ProgramsFiniteMixtures.R 
### use ProgramsAdjacent.R
source("ProgramsFiniteMixtures.R")
source("ProgramsAdjacent.R")


#########################
### Fears data
########################


load("GLES17angst.rda")
dat <- GLES
summary(GLES)
dat <- GLES
summary(GLES)
table(GLES$Unemployment)




####################################
#####  NuclearEnergy and others
#######################################


############## specification 

varloc<-c("Age","Gender","EastWest","Unemployment","Abitur")  ## location variables
varmixt<-c("Age","Gender","EastWest","Unemployment","Abitur")  ## mixture variables
vresp<-"NuclearEnergy"   # response


## alternatives:
#varmixt<-NULL
#vresp<-c("RefugeeCrisis")
#vresp<-c("ClimateChange")

k<-7
deriv<-TRUE

##########################
### without covariate effects on mixture probabilities
##########################

## binomial model
fitbin0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='bin',deriv=TRUE,hessian= TRUE,start=0)
fitbin0
probmixt<-exp(fitbin0$mixture[1,1])/(1+exp(fitbin0$mixture[1,1]))


# with CUB program
CUB<-GEM(Formula(NuclearEnergy~0|Age+Gender++EastWest+Unemployment+Abitur|0), 
            family="cub",data =dat)
summary(CUB)

##  CUBE
CUBE  <-GEM(Formula(NuclearEnergy~0|Age+Gender++EastWest+Unemployment+Abitur|0), 
            family="cube",data =dat)
summary(CUBE)

#not available:
CUBEext  <-GEM(Formula(NuclearEnergy~Age+Gender+Unemployment+EastWest|Age+Gender+Unemployment+EastWest|0), 
            family="cube",data =dat)


## adjacent categories model
fitad0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='adj',deriv=TRUE,hessian= TRUE,start=0)
fitad0
probmixt<-exp(fitad0$mixture[1,1])/(1+exp(fitad0$mixture[1,1]))

## scaled model
fitbinscaled<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='binscaled',deriv,hessian= TRUE,start=0)
fitbinscaled
probmixt<-exp(fitbinscaled$mixture[1,1])/(1+exp(fitbinscaled$mixture[1,1]))

## tests for CUBand binscaled
lr<- -2*(fitbin0$Loglik-fitad0$Loglik)
pval<-1-pchisq(lr, df=k-2)

lr<- -2*(fitbinscaled$Loglik-fitad0$Loglik)
pval<-1-pchisq(lr, df=k-3)

## cumulative model
fitcum0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='cum',deriv,hessian= TRUE,start=0)
fitcum0
probmixt<-exp(fitcum0$mixture[1,1])/(1+exp(fitcum0$mixture[1,1]))



##################
#### with covariates effects on mixture probabilities
########################

varmixt<-varloc
fitad <-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv,hessian= TRUE,start=0,datnew=dat)
fitad 
probmixt<-exp(fitad$mixture[1,1])/(1+exp(fitad$mixture[1,1]))

fitbin <-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv,hessian= TRUE,start=0,datnew=dat)
fitbin 
probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))

fitbinscaled<-Fitmix(vresp,varloc,varmixt,dat,k,type='binscaled',deriv,hessian= TRUE,start=0)
fitbinscaled
probmixt<-exp(fitbinscaled$mixture[1,1])/(1+exp(fitbinscaled$mixture[1,1]))

lr<- -2*(fitbinscaled$Loglik-fitad$Loglik)
pval<-1-pchisq(lr, df=k-3)

### reduced variables

varmixt<-c("Gender","EastWest")
#varmixt<-c("Age")
fitadred <-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv,hessian= TRUE,start=0,datnew=dat)
fitadred 
probmixt<-exp(fitadred$mixture[1,1])/(1+exp(fitadred$mixture[1,1]))

fitbinscaledred<-Fitmix(vresp,varloc,varmixt,dat,k,type='binscaled',deriv,hessian= TRUE,start=0)
fitbinscaledred
probmixt<-exp(fitbinscaledred$mixture[1,1])/(1+exp(fitbinscaledred$mixture[1,1]))


#### dispersion model
fitn<-Fitadj(vresp,varloc,dat,k,deriv='der',hessian,maxit=100,start=0) 
fitn  

vdisp<-varloc
fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian,maxit,start=c(fitn$parameter,rep(0,length(vdisp)))) 
#fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv='der',hessian,maxit,start=0) 
fitnd 

round(fitnd$location, digits=3)
round(fitnd$uncertainty, digits=3)

fitcum <-Fitmix(vresp,varloc,varmixt,dat,k,type='cum',deriv,hessian= TRUE,start=0,datnew=dat)
fitcum

fitcum2 <-Fitmix(vresp,varloc,varmixt,dat,k,type='cum',deriv,hessian= TRUE,start=fitad$parameter,datnew=dat)
fitcum2



#################################
#####  some fits RefugeeCrisis with standard programs
##################################

#### vglm fit
formula<-RefugeeCrisis ~Age+Gender+EastWest+Unemployment

fitvglm <- vglm(formula,family=acat(parallel=TRUE),data=GLES )
fitvglm
summary(fitvglm)

fitvcum <- vglm(formula,family=cumulative(parallel=TRUE),data=GLES )
fitvcum
summary(fitvcum)

### CUB

### without covariates in mixture

CUB0  <-GEM(Formula(RefugeeCrisis~0|Age+Gender+Unemployment+EastWest|0), 
            family="cub",data =GLES)
summary(CUB0)

CUBE  <-GEM(Formula(RefugeeCrisis~0|Age+Gender+Unemployment+EastWest|0), 
            family="cube",data =GLES)
summary(CUBE)



#with covariates

CUB  <-GEM(Formula(RefugeeCrisis~Age+Gender+EastWest+Unemployment|Age+Gender+EastWest+Unemployment|0), 
           family="cub",data =GLES)

summary(CUB) 

CUB0  <-GEM(Formula(NuclearEnergy~0|Age+Gender+Unemployment+EastWest|0), 
            family="cub",data =GLES)
summary(CUB0)

formula<-NuclearEnergy ~Age+Gender+EastWest+Unemployment

fitvglm <- vglm(formula,family=acat(parallel=TRUE),data=GLES )
fitvglm
summary(fitvglm)

fitvcum <- vglm(formula,family=cumulative(parallel=TRUE),data=GLES )
fitvcum
summary(fitvcum)





########################################
#### confidence data:  
###################################


#setwd("C:\\Users\\tutz\\LRZ Sync+Share\\TuRegrMixtureProblems\\R")

load("confidence.rda")

deriv<-TRUE

attach("confidence.rda")
dat<-data
summary(dat)

dat<-dat[dat$income <10000,]  ## outliers
dat$income <-dat$income/1000

d <- density(dat$income) # returns the density data
plot(d)
d <- density(dat$age) # returns the density data
plot(d)

### specification 
varloc<-c("gender","age","interest","income","health")
varmixt<-varloc
vresp<-"justice"
#vresp<-"government"

k<-7

## CUB 

fitbin0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='bin',deriv=TRUE,hessian= TRUE,start=0)
fitbin0
probmixt<-exp(fitbin0$mixture[1,1])/(1+exp(fitbin0$mixture[1,1]))

## CUB with CUB program
CUB  <-GEM(Formula(justice~0|gender+age+interest+income+health|0), 
           family="cub",data =dat)
summary(CUB)

## CUBE
CUBE  <-GEM(Formula(justice~0|gender+age+interest+income+health|0), 
            family="cube",data =dat)
summary(CUBE)

CUBEg  <-GEM(Formula(government~0|gender+age+interest+income+health|0), 
            family="cube",data =dat)
summary(CUBEg)

### adjacent categories
fitad0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='adj',deriv,hessian= TRUE,start=0)
fitad0
probmixt<-exp(fitad0$mixture[1,1])/(1+exp(fitad0$mixture[1,1]))

## cumulative
fitcum0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='cum',deriv,hessian= TRUE,start=0)
fitcum0
probmixt<-exp(fitcum0$mixture[1,1])/(1+exp(fitcum0$mixture[1,1]))




#### with covariates
fitbin <-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv=TRUE,hessian= TRUE,start=0,datnew=dat)
fitbin 
probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))



fitad<-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv=TRUE,hessian=TRUE,start=0)
fitad
probmixt<-exp(fitad$mixture[1,1])/(1+exp(fitad$mixture[1,1]))

fitbinscaled<-Fitmix(vresp,varloc,varmixt,dat,k,type='binscaled',deriv,hessian= TRUE,start=0)
fitbinscaled
probmixt<-exp(fitbinscaled$mixture[1,1])/(1+exp(fitbinscaled$mixture[1,1]))


#### dispersion model
fitn<-Fitadj(vresp,varloc,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=0) 
fitn  

vdisp<-varloc
fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=c(fitn$parameter,rep(0,length(vdisp)))) 
#fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=FALSE,hessian=FALSE,maxit=1,start=0) 
fitnd 

round(fitnd$uncertainty,3)

###other response
vresp<-"healthcare"
healthcare ~ gender+age+interest+income+health

## CUB
CUB  <-GEM(Formula(healthcare~0|gender+age+interest+income+health|0), 
           family="cub",data =dat)
summary(CUB)

CUBE  <-GEM(Formula(healthcare~0|gender+age+interest+income+health|0), 
            family="cube",data =dat)
summary(CUBE)




####################################
###### safety and happiness
#####################################


#setwd("C:\\Users\\tutz\\LRZ Sync+Share\\TuRegAdjCat\\RGeneral")

happy<-readRDS("relgood2M")
summary(happy)
dat<-happy
fitp <- vglm(Safety~age+Gender+Education+Neighbours,
             family=cumulative(parallel=TRUE,reverse=TRUE),
             data=dat)
summary(fitp)

fitp <- vglm(happicat~age+Gender+Education+Neighbours,
             family=cumulative(parallel=TRUE,reverse=TRUE),
             data=dat)
summary(fitp)

varloc<-c("age","Gender","Education","Neighbours")
varmixt<-varloc
vresp<-"Safety"

####happy:
vresp<-"happicat"
k<-10


deriv<-TRUE

fitbin0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='bin',deriv=TRUE,hessian= TRUE,start=0)
fitbin0
probmixt<-exp(fitbin0$mixture[1,1])/(1+exp(fitbin0$mixture[1,1]))


CUBEhappy  <-GEM(Formula(happicat~0|age+Gender+Education+Neighbours|0), 
            family="cube",data =dat)
summary(CUBEhappy)

###  

fitad0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='adj',deriv,hessian= TRUE,start=0)
fitad0
probmixt<-exp(fitad0$mixture[1,1])/(1+exp(fitad0$mixture[1,1]))

fitbinscaled<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='binscaled',deriv,hessian= TRUE,start=0)
fitbinscaled
probmixt<-exp(fitbinscaled$mixture[1,1])/(1+exp(fitbinscaled$mixture[1,1]))



fitcum0<-Fitmix(vresp,varloc,varmixt,dat,k,type='cum',deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.01))
fitcum0 
probmixt<-exp(fitcum0$mixture[1,1])/(1+exp(fitcum0$mixture[1,1]))

CUBE$estimates

## tests
### CUB
lr<- -2*(fitbin0$Loglik-fitad0$Loglik)
pval<-1-pchisq(lr, df=k-2)

#Scaled
lr<- -2*(fitbinscaled$Loglik-fitad0$Loglik)
pval<-1-pchisq(lr, df=k-3)



### with covariates



varmixt<-varloc
fitad<-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv=TRUE,hessian=TRUE,start=0)
fitad
probmixt<-exp(fitad$mixture[1,1])/(1+exp(fitad$mixture[1,1]))

varmixt<-varloc
fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitbin
probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))

fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitbin
probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))

fitmixsc<-Fitmix(vresp,varloc,varmixt,dat,k,type="binscaled",deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitmixsc
probmixt<-exp(fitmixsc$mixture[1,1])/(1+exp(fitmixsc$mixture[1,1]))


#### dispersion model
fitn<-Fitadj(vresp,varloc,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=0) 
fitn  

vdisp<-varloc
fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=c(fitn$parameter,rep(0,length(vdisp)))) 
#fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian,maxit,start=0) 
fitnd 




### happiness  CUB program
CUB  <-GEM(Formula(happicat~0|age+Gender+Education+Neighbours|0), 
           family="cub",data =dat)
summary(CUB)

CUBE  <-GEM(Formula(happicat~0|age+Gender+Education+Neighbours|0), 
            family="cube",data =dat)
summary(CUBE)

## CUB safety
CUB  <-GEM(Formula(Safety~0|age+Gender+Education+Neighbours|0), 
           family="cub",data =dat)
summary(CUB)


CUBESaf  <-GEM(Formula(Safety~0|age+Gender+Education+Neighbours|0), 
               family="cube",data =dat)
summary(CUBESaf)


##########################
### Gender equality (Ingrid)
############################

#etwd("C:\\Users\\tutz\\LRZ Sync+Share\\TuRegrMixtureProblems\\R")
#C:/Users/tutz/LRZ Sync+Share/TuRegrMixtureProblems/R/EVS_Kurzbeschreibung_GT.txt

load("Data.RData")

View(EVS_DE)
dat<-EVS_DE
summary(dat)
dim(dat)

#v72 When a mother works for pay, the children suffer
#v73 A job is alright but what most women really want is a home and children
#v74 All in all, family life suffers when the woman has a full-time job
#v75 A man's job is to earn money; a woman's job is to look after the home and family
#v76 On the whole, men make better political leaders than women do
#v77 A university education is more important for a boy than for a girl
#v78 On the whole, men make better business executives than women do


varloc<-c("male","age","income","townsize")
varmixt<-varloc
vresp<-"v75"


k<-4

#vresp<-"government"

fitbin0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='bin',deriv,hessian= TRUE,start=0)
fitbin0
probmixt<-exp(fitbin0$mixture[1,1])/(1+exp(fitbin0$mixture[1,1]))

fitmixsc0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type="binscaled",deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitmixsc0
probmixt<-exp(fitmixsc0$mixture[1,1])/(1+exp(fitmixsc0$mixture[1,1]))


fitad0<-Fitmix(vresp,varloc,varmixt=NULL,dat,k,type='adj',deriv,hessian= TRUE,start=0)
fitad0
probmixt<-exp(fitad0$mixture[1,1])/(1+exp(fitad0$mixture[1,1]))


################
### ordinal
dat$townsizeord<-dat$townsize
dat$townsizeord<-as.ordered(dat$townsizeord)





varloc<-c("male","age","income","townsizeord")
varmixt<-varloc
vresp<-"v75"



fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitbin
probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))

fitmixsc<-Fitmix(vresp,varloc,varmixt,dat,k,type="binscaled",deriv=TRUE,hessian=TRUE,start=0,lambda=c(0,.0))
fitmixsc
probmixt<-exp(fitmixsc$mixture[1,1])/(1+exp(fitmixsc$mixture[1,1]))
 


fitad<-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv,hessian= TRUE,start=0)
fitad
probmixt<-exp(fitad$mixture[1,1])/(1+exp(fitad$mixture[1,1]))



### for disp generate data
dm<-GenDatSplits(dat,varloc,varmixt)
pred<-dm$pred[,4:7]

dat<-cbind(dat,dm$pred[,4:7])


#### dispersion model
varloc<-c("male","age","income","townsizeord1","townsizeord2","townsizeord3","townsizeord4")
varmixt<-varloc

fitn<-Fitadj(vresp,varloc,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=0) 
fitn  

vdisp<-varloc
fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian=TRUE,maxit=100,start=c(fitn$parameter,rep(0,length(vdisp)))) 
#fitnd<-Fitadjdisp(vresp,varloc,vdisp,dat,k,deriv=TRUE,hessian,maxit,start=0) 
fitnd 













######################################### own programs
#### fits adjacent mixture

#source("ProgramsAdjacent.R")
#source("ProgramsAdjacentCatSpec.R")
#source("ProgramsAdjMixture.R")


varloc<-c("Age","Gender","EastWest","Unemployment")
varmixt<-c("Age","Gender","EastWest","Unemployment")
#varmixt<-NULL
vresp<-c("RefugeeCrisis")
#vresp<-c("ClimateChange")
k<-7



#### with or without start  

fit<-Fitadj(vresp,varloc,dat,k,deriv=TRUE,hessian=TRUE,maxit=500,start=0)
fit

##without start
#parst<-0

##with start
if ( !is.null(varmixt))parst<-c(fit$parameter,1,rep(0.1,length(varmixt)))/2
if ( is.null(varmixt))parst<-c(fit$parameter,3)/2

# fails without penalty, lambda=c(0.0,0.0)
# caution: if lambda not equal c(0.0,0.0) predictors standardized
fitad<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv=TRUE,hessian=TRUE,start=parst,lambda=c(0.01,0.01))
fitad

probmixt<-exp(fitad$mixture[1,1])/(1+exp(fitad$mixture[1,1]))

lr<--2*(fit$Loglik-fitad$Loglik)
if (!is.null(varmixt))pval<-1-pchisq(lr, df=(dim(varmixt)[2]+1))
if ( is.null(varmixt))pval<-1-pchisq(lr, df=1)
lr
pval

round(fit$location, 3)
round(fitad$location, 3)

##########################

##### bin - identical to CUB fit

fitbin<-Fitmixbin(vresp,varloc,varmixt,dat,k,deriv=TRUE,hessian = TRUE,start=0) 
fitbin  



probmixt<-exp(fitbin$mixture[1,1])/(1+exp(fitbin$mixture[1,1]))

round(fitbin$location, 3)
round(fitbin$mixture, 3)
#fitbin<-fitmixbin(resp,k,pred,predmix,deriv='der',hessian = TRUE,maxit,start=0)
#fitbin

fitad$Loglik
fitbin$Loglik

fitad$AIC
fitbin$AIC

### likelihood ratio
lr<--2*(fitbin$Loglik-fitad$Loglik)
pval<-1-pchisq(lr, df=k-2)
lr
pval

#### with covariates in mixture: bin/CUB ok, for adjacent, cumulative less stable 

## CUB
CUB  <-GEM(Formula(RefugeeCrisis~Age+Gender+EastWest+Unemployment|Age+Gender+EastWest+Unemployment|0), 
           family="cub",data =GLES)

summary(CUB) 


#### own 

varloc<-c("Age","Gender","EastWest","Unemployment")
#varmixt<-c("Age")
varmixt<-c("Age","Gender","EastWest","Unemployment")
#varmixt<-c("Age","Gender","EastWest")


fit<-Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit,start)
fit

##without start
#parst<-0

##with start
if ( !is.null(varmixt))parst<-c(fit$parameter,1,rep(0.1,length(varmixt)))/2
if ( is.null(varmixt))parst<-c(fit$parameter,3)/2


fitbin<-Fitmixbin(vresp,varloc,varmixt,dat,k,deriv=TRUE,hessian = TRUE,start=0) 
fitbin  

## alternative, same result
fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv=TRUE,hessian=TRUE,start=0)
fitbin

fitad<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv=TRUE,hessian=TRUE,start=parst,lambda<-c(10.1,0))
fitad

## alternative, same result
fitad<-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv=TRUE,hessian=TRUE,start=0)
fitad

### cumulative:  unemployment bad

fitcum<-Fitcum(vresp,varloc,dat,k,deriv=TRUE,hessian=TRUE,start=0)
fitcum

if ( !is.null(varmixt))parst<-c(fitcum$parameter,1,rep(0.1,length(varmixt)))/2
if ( is.null(varmixt))parst<-c(fitcum$parameter,3)/2

parst<-0
fitcummixt<-Fitmixcum(vresp,varloc,varmixt,dat,k,deriv=TRUE,hessian=TRUE,start=parst)
fitcummixt


############## climate change

CUB  <-GEM(Formula(ClimateChange~Age+Gender+EastWest+Unemployment|Age+Gender+EastWest+Unemployment|0), 
           family="cub",data =GLES)

summary(CUB) 

### own (same for bin)
varloc<-c("Age","Gender","EastWest","Unemployment")
varmixt<-c("Age","Gender","EastWest","Unemployment")

vresp<-c("ClimateChange")

fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv,hessian= TRUE,start=0)
fitbin

###unstable
fitadj<-Fitmix(vresp,varloc,varmixt,dat,k,type='adj',deriv,hessian= TRUE,start=0)
fitadj


#############################
########### climate change


# CUB with covariates

CUB  <-GEM(Formula(ClimateChange~Age+Gender+EastWest+Unemployment|Age+Gender+EastWest+Unemployment|0), 
           family="cub",data =GLES)
summary(CUB) 

CUB0  <-GEM(Formula(ClimateChange~0|Age+Gender+EastWest+Unemployment|0), 
            family="cub",data =GLES)

summary(CUB0)


CUBE1  <-GEM(Formula(ClimateChange~0|Age+Gender+EastWest+Unemployment|0), 
             family="cube",data =GLES)
summary(CUBE1)

### following does not work, CUBE allows only separate variables in the component!
CUBE2  <-GEM(Formula(ClimateChange~Age+Gender+EastWest+Unemployment|
                       Age+Gender+Unemployment+EastWest|0), 
             family="cube",data =GLES)
summary(CUBE2)


### own
hessian<-TRUE
varloc<-c("Age","Gender","EastWest","Unemployment")
varmixt<-c("Age","Gender","EastWest","Unemployment")
#varmixt<-NULL
vresp<-c("ClimateChange")
k<-7

fit<-Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit,start)
fit

##without start
#parst<-0

##with start
if ( !is.null(varmixt))parst<-c(fit$parameter,1,rep(0.1,length(varmixt)))/2
if ( is.null(varmixt))parst<-c(fit$parameter,3)/2

### same as CUB:
fitbin<-Fitmixbin(vresp,varloc,varmixt,dat,k,deriv='der',hessian = TRUE,start=0) 
fitbin  

## alternative, same result
fitbin<-Fitmix(vresp,varloc,varmixt,dat,k,type='bin',deriv,hessian= TRUE,start=0)
fitbin

## deteriorates:
fitad<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv,hessian,start=parst,lambda=c(0,0))
fitad

fitadpen<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv,hessian,start=parst,lambda=c(0.0,0.105))
fitadpen

## cum fits
fitcum<-Fitmixcum(vresp,varloc,varmixt,dat,k,deriv,hessian,start=0)
fitcum

fitcum2<-Fitmixcum(vresp,varloc,varmixt,dat,k,deriv,hessian,start=0,lambda=c(0.0,0.105))
fitcum2


### without covariates

fitad0<-Fitmixadj(vresp,varloc,varmixt=NULL,dat,k,deriv,hessian,start=0)
fitad0
probmixt<-exp(fitad0$mixture[1,1])/(1+exp(fitad0$mixture[1,1]))

fitcum0<-Fitmixcum(vresp,varloc,varmixt=NULL,dat,k,deriv,hessian,start=0)
fitcum0
probmixt<-exp(fitcum0$mixture[1,1])/(1+exp(fitcum0$mixture[1,1]))

Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit=500,start=0)

#### some trials

varloc<-c("Age","Gender","EastWest","Unemployment","Abitur" )
varmixt<-c("Age","Gender","EastWest","Unemployment")
#varmixt<-c("Age")
#varloc<-c("Age")
#varmixt<-NULL
fittr<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv,hessian,start=0)
#fittr<-Fitmixbin(vresp,varloc,varmixt,dat,k,deriv,hessian,start=0)
fittr
exp(fittr$mixture[1,1])/(1+exp(fittr$mixture[1,1]))

