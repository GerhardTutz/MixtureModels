
###  adjacent categories, scaled


#################################
#########  programs basic adjacent categories



Fitadj<-function(vresp,varloc,dat,k,deriv,hessian,maxit,start){  
  
  ### with variables as arguments
  
  # is fitadj<-function(resp,k,pred,disp,deriv,hessian,maxit,start)
  #  but other arguments 
  ### vresp:responses variable
  ### varloc: predictors location 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  
  pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  resp<-as.matrix(as.numeric(dat[,vresp]))
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  
  if(sum(disp)==0){
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv !=TRUE)fitopt <- optim(parst, loglikadjcat, gr = NULL,resp=resp,k=k,pred=pred, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv ==TRUE)fitopt <- optim(parst, loglikadjcat, gr = derloglikadjcat,resp=resp,k=k,pred=pred, method = "BFGS",
                                     lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  #if(sum(disp)!=0){
  #  betadisp<-rep(.01,dim(disp)[2])
  #  parst<-c(beta,beta0,betadisp)
  #  if(sum(start!=0))parst<-start
  #  if (deriv !='der')fitopt <- optim(parst, loglikadjcatdisp, gr = NULL,resp=resp,k=k,pred=pred,disp=disp, method = "Nelder-Mead",
  #                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  #  if (deriv =='der')fitopt <- optim(parst, loglikadjcatdisp, gr = derloglikadjcatdisp,resp=resp,k=k,pred=pred,disp=disp, method = "BFGS",
  #                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  #  }
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  pardisp<-0
  if(sum(disp!=0)) pardisp<-fitopt$par[(p+k):(p+k+dim(disp)[2]-1)]
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){
    #eps <- 1e-6
    #Hreg <- fitopt$hessian + diag(eps, nrow(fitopt$hessian))
    #hessinv<-ginv(Hreg)  
  hessinv<-ginv(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  zval<-fitopt$par[1:p]/stderr
  pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
  
  if(sum(disp!=0)){stderrdisp<-std[(p+k):(p+k+dim(disp)[2]-1)]
  zvaldisp<-pardisp/stderrdisp
  pardisp<-cbind(pardisp,stderrdisp,zvaldisp)}
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], stderr,zval,pval)
  
  location<-as.data.frame(location)
  
  names(location)[1] <- "Estimates"
  names(location)[2] <- "std err"
  names(location)[3] <- "z-values"
  names(location)[4] <- "p-values"
  row.names(location)<-varloc 
  
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "parunc"= pardisp,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

##############################################



Fitadjdisp<-function(vresp,varloc,vdisp,dat,k,deriv,hessian,maxit,start){  
  
  ### with variables as arguments
  
  # is fitadj<-function(resp,k,pred,disp,deriv,hessian,maxit,start)
  #  but other arguments 
  ### vresp:responses variable
  ### varloc: predictors location 
  ### disp: dispersion variables
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  
  
  pred <-as.matrix(dat[,varloc])
  disp <-as.matrix(dat[,vdisp])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  resp<-as.matrix(as.numeric(dat[,vresp]))
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  
  if(is.null(vdisp )){
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv !=TRUE)fitopt <- optim(parst, loglikadjcat, gr = NULL,resp=resp,k=k,pred=pred, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv ==TRUE)fitopt <- optim(parst, loglikadjcat, gr = derloglikadjcat,resp=resp,k=k,pred=pred, method = "BFGS",
                                     lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  if(!is.null(vdisp)){
    betadisp<-rep(.01,length(vdisp))
    parst<-c(beta,beta0,betadisp)
    if(sum(start!=0))parst<-start
    if (deriv !=TRUE)fitopt <- optim(parst, loglikadjcatdisp, gr = NULL,resp=resp,k=k,pred=pred,disp=disp, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv ==TRUE)fitopt <- optim(parst, loglikadjcatdisp, gr = derloglikadjcatdisp,resp=resp,k=k,pred=pred,disp=disp, method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  uncertainty<-NULL
  
  if(hessian ==TRUE){hessinv<-ginv(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  #stderr<-std[1:p]
  zval<-fitopt$par/std
  pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
  
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], std[1:p],zval[1:p],pval[1:p])
  
  location<-as.data.frame(location)
  
  names(location)[1] <- "Estimates"
  names(location)[2] <- "std err"
  names(location)[3] <- "z-values"
  names(location)[4] <- "p-values"
  row.names(location)<-varloc 
  
  if(is.null(vdisp))uncertainty<-NULL
  if(!is.null(vdisp)) {pardisp<-fitopt$par[(p+k):(p+k+length(vdisp)-1)]
  uncertainty<-cbind( pardisp, std[(p+k):(p+k+length(vdisp)-1)],zval[(p+k):(p+k+length(vdisp)-1)],pval[(p+k):(p+k+length(vdisp)-1)])
  uncertainty<-as.data.frame(uncertainty)
  
  names(uncertainty)[1] <- "Estimates"
  names(uncertainty)[2] <- "std err"
  names(uncertainty)[3] <- "z-values"
  names(uncertainty)[4] <- "p-values"
  row.names(uncertainty)<-vdisp 
  
  }
  
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "uncertainty"= uncertainty,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

##############################################












###############################
fitadj<-function(resp,k,pred,disp,deriv,hessian,maxit,start){
  
  ### with matrices as arguments
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### disp: predictors in dispersion part for uncertainty model
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### maxit: number of iterations
  ### start: start values, none if 0
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  
  if(sum(disp)==0){
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv ==FALSE)fitopt <- optim(parst, loglikadjcat, gr = NULL,resp=resp,k=k,pred=pred, method = "Nelder-Mead",
                  lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)

  if (deriv ==TRUE)fitopt <- optim(parst, loglikadjcat, gr = derloglikadjcat,resp=resp,k=k,pred=pred, method = "BFGS",
                  lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  }
  
  if(sum(disp!=0)){
    betadisp<-rep(.01,dim(disp)[2])
    parst<-c(beta,beta0,betadisp)
    if(sum(start!=0))parst<-start
    if (deriv !='der')fitopt <- optim(parst, loglikadjcatdisp, gr = NULL,resp=resp,k=k,pred=pred,disp=disp, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt <- optim(parst, loglikadjcatdisp, gr = derloglikadjcatdisp,resp=resp,k=k,pred=pred,disp=disp, method = "BFGS",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  pardisp<-0
  if(sum(disp!=0)) pardisp<-fitopt$par[(p+k):(p+k+dim(disp)[2]-1)]
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){hessinv<-ginv(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std[1:p]
    zval<-fitopt$par[1:p]/stderr
        if(sum(disp!=0)){stderrdisp<-std[(p+k):(p+k+dim(disp)[2]-1)]
        zvaldisp<-pardisp/stderrdisp
        pardisp<-cbind(pardisp,stderrdisp,zvaldisp)}
      
    #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
    location<-cbind( fitopt$par[1:p], stderr,zval)
            }
  
  
  
   #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "parunc"= pardisp,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

  
##################################
loglikadjcat<-function(par,resp,k,pred){

### negative loglikelihood  
 # par is betavector,beta_02,... 

#zgamma<- 1  # not yet needed'  
    
p<-dim(pred)[2]
beta<-par[1:p]  ### first are covariate weights
beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds

###### probabilities

n<-dim(resp)[1]

prob<- matrix(0,n,k)



loglik<-0
for(i in 1:n){

### etaterm
  eta<- matrix(0,1,k) 
  for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)}

### etasumterm
  etasum<- matrix(0,1,k)
  for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  

### expterm  
expterm<- matrix(0,1,k)
for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}

#  sum denominator
sum <-1
for (r in 2:k){sum<-sum+expterm[1,r]}

for (r in 1:k){prob[i,r]<-expterm[1,r]/sum}
sum(prob[i,])
loglik<-loglik+log(prob[i,])[resp[i]]
#print(i)
#print(log(prob[i,]))
}
#log(prob)

##### log-lik

loglikneg<- -loglik

newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
return(newList)}



##################################
derloglikadjcat<-function(par,resp,k,pred){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  der2<-matrix(0,(k-1),1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-der1+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der[1:p,1]<-der[1:p,1]+(cat-1)*(as.matrix(pred[i,]))-dersum[,1]/sum
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
  }  ## end i
  
  ##### log-lik
  der<-c(der1,der2)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

###############################################



loglikadjcatdisp<-function(par,resp,k,pred,disp){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # disp is matrix with dispersion covariates
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-(beta0[r]+(pred[i,]%*%beta))/exp(disp[i,]%*%betadisp)}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum}
    sum(prob[i,])
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}


##################################
derloglikadjcatdisp<-function(par,resp,k,pred,disp){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  der2<-matrix(0,(k-1),1)
  der3<-matrix(0,pdisp,1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-(beta0[r]+(pred[i,]%*%beta))/exp(disp[i,]%*%betadisp)}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    sum1<-0
    sum3<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) {sum1<-sum1 +prob[i,r]*(r-1) 
              sum3<-sum3+ prob[i,r]*etasum[1,r]    }
            }
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-der1+((cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,])))*as.numeric(exp(-disp[i,]%*%betadisp))
    #der[1:p,1]<-der[1:p,1]+(cat-1)*(as.matrix(pred[i,]))-dersum[,1]/sum
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-exp(-disp[i,]%*%betadisp)
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])*exp(-disp[i,]%*%betadisp)
    }
   # dispersion
    der3<-der3+etasum[1,cat]*(-as.matrix(disp[i,])) +sum3*as.matrix(disp[i,])
    
    
    }  ## end i
  
  ##### log-lik
  der<-c(der1,der2,der3)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

###############################################




loglikadjcatRepeat<-function(par,resp,k,pred){
  
  ### negative loglikelihood  Repeat - maybe redundant
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum}
    sum(prob[i,])
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

###############################################

