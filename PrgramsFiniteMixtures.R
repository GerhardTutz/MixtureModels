

###########

##### Programs Finite Mixtures  

################################

#### ordered predictors automatically transformed into split variables
#### group lasso included
#### fitted:   adjacent categories, binomial, binomial scaled, cumulative mixture models


### wrapper Fitmix calls binomial, adjacent categories, cumulative programs

###ProgramsAdjacent needed

library("VGAM")
library("xtable")
library("CUB")
install.packages("FastCUB")
library("FastCUB")
library(MASS)


Fitmix<-function(vresp,varloc,varmixt,dat,k,type,deriv=NULL,hessian,start,lambda =NULL,datnew=NULL, stand=NULL){
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  ### lambda: vector penalty parameter localisation, mixture, e.g. lambda<-c(.1,1)
  ## type: 'bin' for binomial
  #        'binscaled' for binomial scaled
  #        'adj' for adjacent categories
  #        'cum' for cumulative
  
  if(is.null(deriv))deriv <-'der'
  
  if(type =='bin')fit<-Fitmixbin(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda,datnew, stand)
  if(type =='adj')fit<-Fitmixadj(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda,datnew, stand)
  if(type =='cum'){fit<-
    Fitmixcum(vresp = vresp,varloc = varloc,varmixt = varmixt,dat = dat,k = k,deriv = deriv,
              hessian = hessian,start = start,lambda = lambda,datnew, stand)}
  
  if(type =='binscaled')fit<-FitmixbinScaled(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda,datnew, stand)
  return(fit)
}



###############################
#fitmixadj<-function(resp,k,pred,predmix,deriv,hessian,maxit,start){

Fitmixadj<-function(vresp,varloc,varmixt,dat,k,deriv=NULL,hessian,start,lambda=NULL,datnew=NULL, stand=NULL){
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  ### standardize: if TRUe predictors standardized
  ###  lassoparameters, if non equql (0,0) predictors standardized
  
  ### needs fitadj
  disp<-0  ## obsolete value
  if(is.null(lambda ))lambda <-c(0,0)
  if(is.null(deriv ))deriv <-'der'
  
  #pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  #resp<-as.matrix(as.numeric(dat[,vresp]))
  
  ### design for ordered variables
  
  #all(!sapply(dat, is.ordered)) ## TRUE if all are false
  
  dm<-GenDatSplits(dat,varloc,varmixt)
  pred<-dm$pred
  predmix<-dm$predmix
  namespred <-dm$namespred 
  namespredmix<-dm$namespredmix
  dm$locvar  ## list with placements for variables
  dm$mixvar
  listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                   mixvar=dm$namespredmix, plmixvar=dm$mixvar)
  #listlambda$locvar[[1]]
  
  resp<-as.matrix(as.numeric(dat[,vresp]))
  
  maxit<-500
  
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  
  
  if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
  if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
  
  ###### standardization
  if (!is.null(stand)){
    pred<-scale(pred)
    if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
  }
  ####################
  
  pmix<-dim(predmix)[2]
  
  gama<-NULL
  if(dim(predmix)[2]>1)gama<-c(0,rep(.01,dim(predmix)[2]-1))
  
  ##start values
  
  if(sum(start)==0){
    fit<-fitadj(resp,k,pred,disp=0,deriv=TRUE,hessian,maxit,start=0)
    gama<-c(0.5,rep(.01,pmix-1))
    parst<-c(fit$parameter,gama)
  }
  if(sum(start!=0))parst<-start
  
  ## fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtadjcat, gr = NULL,resp=resp,k=k,pred=pred,
                                    predmix=predmix,listlambda=listlambda, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                                    hessian = hessian)
  
  if (deriv ==TRUE)fitopt <- optim(parst, loglikmixtadjcat, gr = derloglikmixtadjcat,resp=resp,k=k,
                                   pred=pred,predmix=predmix,listlambda=listlambda, method = "BFGS",
                                   lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                                   hessian = hessian)
  
  listn<-listlambda
  listn[[1]]<-c(0,0)
  loglikunpen<--loglikmixtadjcat(fitopt$par,resp,k,pred,predmix,listn)
  loglikpen<--fitopt$value
  der<-derloglikmixtadjcat(fitopt$par,resp,k,pred,predmix,listlambda)
  #loglikmixtadjcat(parst,resp,k,pred,predmix,listn)
  #penlass(.001,parst,rep(1,length(par)))
  #derlass(.001,parst,rep(1,length(par)))
  
  AIC<- -2*(loglikunpen-length(parst))
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    eps <- 1e-6
    Hreg <- fitopt$hessian + diag(eps, nrow(fitopt$hessian))
    hessinv<-ginv(Hreg)
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    
    location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p],pval[1:p])
    mixture<-cbind(fitopt$par[(p+k):length(parst)], stderr[(p+k):length(parst)],
                   zval[(p+k):length(parst)],pval[(p+k):length(parst)])
    
    
    location<-as.data.frame(location)
    
    names(location)[1] <- "Estimates"
    names(location)[2] <- "std err"
    names(location)[3] <- "z-values"
    names(location)[4] <- "p-values"
    row.names(location)<-namespred #c(seq(1,(knew-1),1),paste(pred))
    
    mixture<-as.data.frame(mixture)
    
    names(mixture)[1] <- "Estimates"
    names(mixture)[2] <- "std err"
    names(mixture)[3] <- "z-values"
    names(mixture)[4] <- "p-values"
    row.names(mixture)<-c("const",namespredmix)
  }
  
  if(hessian ==FALSE){
    location<- fitopt$par[1:p]
    mixture<-fitopt$par[(p+k):length(parst)]
  }
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  
  #### prediction
  loglikpred<-0
  
  if(!is.null(datnew )){
    dm<-GenDatSplits(datnew,varloc,varmixt)
    pred<-dm$pred
    predmix<-dm$predmix
    namespred <-dm$namespred 
    namespredmix<-dm$namespredmix
    dm$locvar  ## list with placements for variables
    dm$mixvar
    listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                     mixvar=dm$namespredmix, plmixvar=dm$mixvar)
    
    resp<-as.matrix(as.numeric(datnew[,vresp]))
    
    if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
    if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
    
    ###### standardization
    if (sum(lambda[[1]])>0){
      pred<-scale(pred)
      if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
    }
    ####################
    listn<-listlambda
    listn[[1]]<-c(0,0)
    
    loglikpred<--loglikmixtadjcat(fitopt$par,resp,k,pred,predmix,listn) 
    
  }
  
  newList <- list("Loglik"= loglikunpen,"LoglikPen"= loglikpen,"AIC"=AIC, "location"=location, "mixture"= mixture,
                  "parameter"=fitopt$par, "convergence"= fitopt$convergence,
                  "der"=der,"loglikpred"=loglikpred)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

######################################
GenDatSplits<-function(dat,varloc,varmixt){
  
  ### generates data as matrix and names neede if ordered factors are present
  ### generates list of placements for variables
  
  loc<-makeDesign0(dat,varloc)
  loc$X <- loc$X[, colnames(loc$X) != "(Intercept)", drop = FALSE]
  pred <-as.matrix(loc$X)
  namespred<-names(loc$X)
  #loc$variable_index_list[[1]]   ### list with categories
  
  if ( is.null(varmixt)){
    predmix<-NULL
    namespredmix<-NULL
  } else {
    mix<-makeDesign0(dat,varmixt)
    mix$X <- mix$X[, colnames(mix$X) != "(Intercept)", drop = FALSE]
    predmix <-as.matrix(mix$X)
    namespredmix<-names(mix$X)
    #mix$variable_index_list[[2]]  ### list with categories
  }
  
  if ( !is.null(varmixt)){newList <- list("pred"= pred,"predmix"= predmix,"namespred"=namespred,"namespredmix"=namespredmix,"locvar"=loc$variable_index_list,
                  "mixvar"=mix$variable_index_list)
  } else {  
   newList <- list("pred"= pred,"predmix"= predmix,"namespred"=namespred,"namespredmix"=namespredmix,"locvar"=loc$variable_index_list)
   }  

return(newList)
  }
######################
 

###########################

makeDesign <- function(dat, var){
  ### alte Version!
  ord_vars <- sapply(dat, is.ordered)
  dat[ord_vars] <- lapply(dat[ord_vars], function(x) {
    contrasts(x) <- contr.split(nlevels(x))
    x
  })
  
  all <- reformulate(var)   # keep intercept
  X <- model.matrix(all, data = dat)
  
  newList <- list("X" = as.data.frame(X),
                  "names" = colnames(X))
  
  return(newList)
}

############################################

contr.split <- function(n) {
  mat <- matrix(0, n, n - 1)
  for (i in 2:n) {
    mat[i:n, i - 1] <- 1
  }
  mat
}
#########################################
##############################
makeDesign0 <- function(dat, var){
  
  ord_vars <- sapply(dat, is.ordered)
  
  dat[ord_vars] <- lapply(dat[ord_vars], function(x) {
    contrasts(x) <- contr.split(nlevels(x))
    x
  })
  
  f <- reformulate(var)
  X <- model.matrix(f, data = dat)
  
  assign_vec <- attr(X, "assign")
  
  # remove intercept columns completely
  keep <- assign_vec != 0
  X <- X[, keep, drop = FALSE]
  assign_vec <- assign_vec[keep]
  
  # term labels in correct order
  term_labels <- attr(terms(f), "term.labels")
  
  variable_index_list <- lapply(seq_along(term_labels), function(i) {
    which(assign_vec == i)
  })
  
  names(variable_index_list) <- var
  
  list(
    X = as.data.frame(X),
    names = colnames(X),
    variable_index_list = variable_index_list
  )
}
#############

##################################
loglikmixtadjcat<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  #par<-parst
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds  with o
  gama<-par[(p+k):(p+k+dim(predmix)[2]-1)]
  
  pred <- as.matrix(pred)
  predmix <- as.matrix(predmix)
  ###### probabilities
  
  n<-dim(pred)[1]
  
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){prob[i,r]<-pimix*expterm[1,r]/sumn+ (1-pimix)/k}   
    #sum(prob[i,])
    #sum(expterm[i,])
    
    
    
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ### pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    #varn<-length(listlambda$locvar)   
    varn<-length(listlambda$pllocvar)
    if (varn >0){
      for (v in 1:varn)  
        loglik<-loglik-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*NormShiftv(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        loglik<-loglik-sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*NormShiftv(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  ##### log-lik
  
  loglikneg<- -loglik
  
  #newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(loglikneg)}

##################################

###################
derloglikmixtadjcat<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  #print(par)
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  gama<-par[(p+k):(p+k+dim(predmix)[2]-1)]
  
  pred <- as.matrix(pred)
  predmix <- as.matrix(predmix)
  
  n<-dim(pred)[1]
  
  prob<- matrix(0,n,k)
  mixtprob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  der2<-matrix(0,(k-1),1)
  der1mixt<-matrix(0,p+(k-1),1)
  der2mixt<-matrix(0,dim(predmix)[2],1)
  
  der<-matrix(0,p+(k-1)+dim(predmix)[2],1)
  #dim(der)
  #loglik<-0
  
  for(i in 1:n){
    
    ###### probabilities adj categories
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sumn 
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    sum(prob[i,])
    
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der[1:p,1]<-der[1:p,1]+(cat-1)*(as.matrix(pred[i,]))-dersum[,1]/sum
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-ind1-sum(prob[i,j:k])
    }
    deradj<-c(der1,der2)
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){mixtprob[i,r]<-pimix*prob[i,r]+ (1-pimix)/k}
    mixtsel<-mixtprob[i,][resp[i]]
    #sum(prob[i,])
    #sum(mixtprob[i,])
    der1mixt<-der1mixt+as.numeric(pimix*prob[i,][resp[i]])*matrix(deradj)/mixtsel
    
    der2mixt<-der2mixt+as.numeric((prob[i,][resp[i]]-1/k)*pimix*(1-pimix))*matrix(predmix[i,])/
      (mixtsel)
    
  }## end i
  
  
  
  
  ##### der
  der<-c(der1mixt,der2mixt)
  #par[1:p]
  
  ### pen
  #old
  #if (lambda[1] >0)der[1:p]<-der[1:p]-derlass(lambda[1],par[1:p],rep(1,length(beta)))
  #if ((lambda[2] >0)&(length(gama)>1))der[(p+k+1):(p+k+dim(predmix)[2]-1)]<-
  #  der[(p+k+1):(p+k+dim(predmix)[2]-1)]-derlass(lambda[2],gama[2:length(gama)],rep(1,length(gama)-1))
  
  ### new pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        der[listlambda$pllocvar[[v]]]<-der[listlambda$pllocvar[[v]]]-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*
          NormShiftvder(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        der[p+k-1+listlambda$plmixvar[[v]]+1]<-der[p+k-1+listlambda$plmixvar[[v]]+1]-
          sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*
          NormShiftvder(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

###############################################

derloglikmixtadjcatchat <- function(par, resp, k, pred, predmix){
  
  p <- dim(pred)[2]
  beta <- par[1:p]
  beta0 <- c(0, par[(p+1):(p+k-1)])
  gama <- par[(p+k):(p+k+dim(predmix)[2]-1)]
  
  n <- dim(pred)[1]
  
  der_beta  <- rep(0, p)
  der_beta0 <- rep(0, k-1)
  der_gamma <- rep(0, dim(predmix)[2])
  
  for(i in 1:n){
    
    # ---- linear predictor (scalar!) ----
    linpred <- as.numeric(pred[i,] %*% beta)
    
    # ---- eta ----
    eta <- numeric(k)
    for (r in 2:k){
      eta[r] <- beta0[r] + linpred
    }
    
    # ---- cumulative sums ----
    etasum <- numeric(k)
    for (r in 2:k){
      etasum[r] <- sum(eta[2:r])
    }
    
    # ---- probabilities (adjacent categories) ----
    expterm <- exp(etasum)
    sumn <- sum(expterm)
    p_adj <- expterm / sumn
    
    # ---- mixture ----
    pimix <- as.numeric(plogis(predmix[i,] %*% gama))
    p_mix <- pimix * p_adj + (1 - pimix)/k
    
    r_obs <- resp[i]
    P_i <- p_mix[r_obs]
    
    # ---- expectation term ----
    E_r_minus1 <- sum(p_adj * (0:(k-1)))
    
    # ---- score for beta ----
    score_beta <- (r_obs - 1 - E_r_minus1) * pred[i,]
    
    # ---- score for thresholds ----
    score_beta0 <- numeric(k-1)
    for(j in 2:k){
      ind <- as.numeric(r_obs >= j)
      score_beta0[j-1] <- ind - sum(p_adj[j:k])
    }
    
    # ---- mixture weight ----
    w <- pimix / P_i
    
    der_beta  <- der_beta  + w * score_beta
    der_beta0 <- der_beta0 + w * score_beta0
    
    # ---- gamma derivative ----
    der_gamma <- der_gamma +
      ((p_adj[r_obs] - 1/k) / P_i) *
      pimix * (1 - pimix) * predmix[i,]
  }
  
  # ---- combine + NEGATE (because we minimize) ----
  der <- -c(der_beta, der_beta0, der_gamma)
  
  return(der)
}

###################claude
derloglikmixtadjcatCl <- function(par, resp, k, pred, predmix){
  
  p <- dim(pred)[2]
  beta  <- par[1:p]
  beta0 <- c(0, par[(p+1):(p+k-1)])
  gama  <- par[(p+k):(p+k+dim(predmix)[2]-1)]
  
  n <- dim(pred)[1]
  
  prob     <- matrix(0, n, k)
  mixtprob <- matrix(0, n, k)
  
  der1mixt <- matrix(0, p+(k-1), 1)
  der2mixt <- matrix(0, dim(predmix)[2], 1)
  
  for(i in 1:n){
    
    # eta terms
    eta <- matrix(0, 1, k)
    for (r in 2:k){ eta[1,r] <- beta0[r] + (pred[i,] %*% beta) }
    
    # cumulative eta sums
    etasum <- matrix(0, 1, k)
    for (r in 2:k){ etasum[1,r] <- sum(eta[1,2:r]) }
    
    # exp terms
    expterm <- matrix(0, 1, k)
    for (r in 1:k){ expterm[1,r] <- exp(etasum[1,r]) }
    
    # denominator — avoid overwriting base sum()
    sumn <- 1
    for (r in 2:k){ sumn <- sumn + expterm[1,r] }
    
    # adj-category probabilities and E[r-1] = sum1
    sum1 <- 0
    for (r in 1:k){
      prob[i,r] <- expterm[1,r] / sumn
      if(r > 1) sum1 <- sum1 + prob[i,r] * (r-1)
    }
    
    cat <- resp[i]
    
    # gradient w.r.t. beta (covariate weights)
    der1 <- (cat-1) * as.matrix(pred[i,]) - sum1 * as.matrix(pred[i,])
    
    # gradient w.r.t. beta0 (thresholds)
    der2 <- matrix(0, k-1, 1)
    for (j in 2:k){
      ind1 <- if(cat >= j) 1 else 0
      der2[j-1] <- ind1 - sum(prob[i, j:k])
    }
    
    deradj <- rbind(as.matrix(der1), as.matrix(der2))  # (p+k-1) x 1
    
    # mixture probability
    pimix <- exp(predmix[i,] %*% gama) / (1 + exp(predmix[i,] %*% gama))
    
    for (r in 1:k){ mixtprob[i,r] <- pimix * prob[i,r] + (1-pimix)/k }
    
    mixtsel <- mixtprob[i, cat]   # scalar: P_mix for observed category
    
    # --- FIX: correct score factor is pimix/mixtsel, not pimix*mixtsel/mixtsel ---
    der1mixt <- der1mixt + as.numeric(pimix / mixtsel) * deradj
    
    # gradient w.r.t. gama (mixing weights) — was already correct
    der2mixt <- der2mixt +
      as.numeric((prob[i, cat] - 1/k) * pimix * (1-pimix)) *
      as.matrix(predmix[i,]) / mixtsel
    
  }  # end i
  
  der <- -c(der1mixt, der2mixt)   # negative gradient (for minimization)
  return(der)
}

###############################
fitmixbin<-function(resp,k,pred,predmix,deriv,hessian,maxit,start,lambda=NULL){
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### maxit: number of iterations
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  maxit<-500
  
  pmix<-dim(predmix)[2]
  
  ##start values
  
  if(sum(start)!=0){
    fit<-fitadj(resp,k,pred,disp=0,deriv='der',hessian= FALSE,maxit,start=0)
    fit
    if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
    if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
    gama<-NULL
    if(dim(predmix)[2]>1)gama<-c(0,rep(.01,dim(predmix)[2]-1))
    parst<-c(fit$parameter[1:p]/2,0,gama)
  }
  if(sum(start)!=0)parst<-start
  
  ## fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv ==TRUE){
    fitoptnd <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix, method = "Nelder-Mead",
                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    fitopt <- optim(fitoptnd$par, loglikmixtbin, gr = derloglikmixtbin,resp=resp,k=k,
                    pred=pred,predmix=predmix, method = "BFGS",
                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                    hessian = hessian)
  }
  
  #loglikmixtadjcat(parst,resp,k,pred,predmix)
  #derloglikmixtbin(fitopt$par,resp,k,pred,predmix)
  #derloglikmixtbin(parst,resp,k,pred,predmix)
  
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    hessinv<-ginv(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    
    location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p])
    mixture<-cbind(fitopt$par[(p+2):length(parst)], stderr[(p+2):length(parst)],zval[(p+2):length(parst)])
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location, "mixture"= mixture,
                  "parunc"= pardisp,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}


#######################################################

loglikmixtbin<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  
  ### thresholds
  thr<-NULL
  for (r in 2:k) thr<-c(thr,log((k-r+1)/(r-1)))
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds  with o
  gama<-par[(p+2):length(par)]
  
  #dim(predmix)
  #length(par)
  ###### probabilities
  
  n<-dim(pred)[1]
  
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){prob[i,r]<-pimix*expterm[1,r]/sumn+ (1-pimix)/k}
    #sum(prob[i,])
    
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ### pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        loglik<-loglik-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*NormShiftv(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        loglik<-loglik-sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*NormShiftv(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  #newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(loglikneg)}

##################################

###################
derloglikmixtbin<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  
  ### thresholds
  thr<-NULL
  for (r in 2:k) thr<-c(thr,log((k-r+1)/(r-1)))
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds   
  gama<-par[(p+2):length(par)]
  
  
  n<-dim(pred)[1]
  
  prob<- matrix(0,n,k)
  mixtprob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  #der2<-matrix(0,1,1)
  der1mixt<-matrix(0,p+1,1)  # for beta, beta_0
  der2mixt<-matrix(0,dim(predmix)[2],1) # for gamma
  
  der<-matrix(0,p+1+dim(predmix)[2],1)
  #dim #
  #loglik<-0
  
  for(i in 1:n){
    
    ###### probabilities adj categories
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sumn
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    sum(prob[i,])
    
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der1 <- as.matrix(((cat-1) - sum1) * pred[i,])
    
    ### thresholds
    
    #for (j in 2:k){
    #  ind1<-0
    #  if(cat>=j)ind1<-1
    #  der2[j-1]<-ind1-sum(prob[i,j:k])    }
    
    #sumdum<-0
    #for (j in 2:k)sumdum<-sumdum+prob[i,j]*(j-1)
    #der2<-cat-1-sumdum
    
    der2 <- 0
    for (j in 2:k){
      ind1 <- 0
      if(cat >= j) ind1 <- 1
      der2 <- der2 + ind1 - sum(prob[i, j:k])
    }
    
    
    deradj<-c(der1,der2)
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){mixtprob[i,r]<-pimix*prob[i,r]+ (1-pimix)/k}
    mixtsel<-mixtprob[i,][resp[i]]
    #sum(prob[i,])
    #sum(mixtprob[i,])
    der1mixt<-der1mixt+as.numeric(pimix*prob[i,][resp[i]])*matrix(deradj)/mixtsel
    
    der2mixt<-der2mixt+as.numeric((prob[i,][resp[i]]-1/k)*pimix*(1-pimix))*matrix(predmix[i,])/
      (mixtsel)
    
  }## end i
  
  #####  
  der<-c(der1mixt,der2mixt)
  
  ### new pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        der[listlambda$pllocvar[[v]]]<-der[listlambda$pllocvar[[v]]]-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*
          NormShiftvder(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        der[p+1+listlambda$plmixvar[[v]]+1]<-der[p+1+listlambda$plmixvar[[v]]+1]-
          sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*
          NormShiftvder(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  
  
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}
#############################

Fitmixbin <-function(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda=NULL,datnew=NULL, stand=NULL){
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  if(is.null(lambda ))lambda <-c(0,0)
  #pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  ### design for ordered variables
  
  dm<-GenDatSplits(dat,varloc,varmixt)
  pred<-dm$pred
  predmix<-dm$predmix
  namespred <-dm$namespred 
  namespredmix<-dm$namespredmix
  dm$locvar  ## list with placements for variables
  dm$mixvar
  listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                   mixvar=dm$namespredmix, plmixvar=dm$mixvar)
  #listlambda$locvar[[1]]
  
  
  
  resp<-as.matrix(as.numeric(dat[,vresp]))
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
  if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
  
  ###### standardization
  if (!is.null(stand)){
    pred<-scale(pred)
    if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
  }
  ####################
  maxit<-500
  
  pmix<-dim(predmix)[2]
  
  ##start values
  
  #fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = FALSE,maxit=maxit,start=0)
  #fit
  if(sum(start)==0){
    #fit<-Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit,start=0)
    gama<-NULL
    if(dim(predmix)[2]>1)gama<-c(0.5,rep(.01,dim(predmix)[2]-1))
    if(dim(predmix)[2]==1)gama<-0.5
    parst<-c(rep(0,p),0,gama) ###!
  }
  
  if(sum(start)!=0)parst<-start
  
  ## fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix,
                                    listlambda=listlambda,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv ==TRUE){
    if(sum(start)==0){maxit<-100
    fitoptnd <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix, 
                      listlambda=listlambda,method = "Nelder-Mead",
                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    parst<-fitoptnd$par}
    
    fitopt <- optim(parst, loglikmixtbin, gr = derloglikmixtbin,resp=resp,k=k,
                    pred=pred,predmix=predmix,listlambda=listlambda, method = "BFGS",
                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                    hessian = hessian)
  }
  
  #loglikmixtadjcat(parst,resp,k,pred,predmix)
  der<-derloglikmixtbin(fitopt$par,resp,k,pred,predmix,listlambda)
  #derloglikmixtbin(parst,resp,k,pred,predmix)
  listn<-listlambda
  listn[[1]]<-c(0,0)
  loglikunpen<--loglikmixtbin(fitopt$par,resp,k,pred,predmix,listn)
  loglikpen<--fitopt$value
  
  
  AIC<- -2*(loglikunpen-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  mixture<-fitopt$par[(p+2):length(parst)]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    #hessinv <-  ginv(fitopt$hessian)
    H <- fitopt$hessian
    eps <- 1e-8
    
    hessinv<- solve(H + diag(eps, nrow(H))) 
    
    #hessinv<-ginv(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p],pval[1:p])
    mixture<-cbind(fitopt$par[(p+2):length(parst)], stderr[(p+2):length(parst)],
                   zval[(p+2):length(parst)],pval[(p+2):length(parst)])
    
    location<-as.data.frame(location)
    
    names(location)[1] <- "Estimates"
    names(location)[2] <- "std err"
    names(location)[3] <- "z-values"
    names(location)[4] <- "p-values"
    row.names(location)<-namespred 
    
    mixture<-as.data.frame(mixture)
    
    names(mixture)[1] <- "Estimates"
    names(mixture)[2] <- "std err"
    names(mixture)[3] <- "z-values"
    names(mixture)[4] <- "p-values"
    row.names(mixture)<-c("const",namespredmix)
    
  }
  
  #### prediction
  loglikpred<-0
  
  if(!is.null(datnew )){
    dm<-GenDatSplits(datnew,varloc,varmixt)
    pred<-dm$pred
    predmix<-dm$predmix
    namespred <-dm$namespred 
    namespredmix<-dm$namespredmix
    dm$locvar  ## list with placements for variables
    dm$mixvar
    listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                     mixvar=dm$namespredmix, plmixvar=dm$mixvar)
    
    resp<-as.matrix(as.numeric(datnew[,vresp]))
    
    if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
    if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
    
    ###### standardization
    if (sum(lambda[[1]])>0){
      pred<-scale(pred)
      if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
    }
    ####################
    listn<-listlambda
    listn[[1]]<-c(0,0)
    
    loglikpred<- -loglikmixtbin(fitopt$par,resp,k,pred,predmix,listn) 
    
  }
  
  
  newList <- list("Loglik"= loglikunpen,"LoglikPen"= loglikpen,"AIC"=AIC, "location"=location, "mixture"= mixture,
                  "parameter"=fitopt$par, "convergence"= fitopt$convergence,"der"=der,"loglikpred"=loglikpred)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

#################################


#############################
penlass<-function(lambdal,par,weights){
  pen<-0
  c<-.001
  #for (v in glob) pen<- pen +NormShiftv(beta[v],c)
  for (v in 1:length(par)) pen<- pen +weights[v]*NormShiftv(par[v],c)
  #loglik<-loglik- lambdagl*pen
  penterm<-lambdal*pen
  return(penterm)
}
###########################
derlass<-function(lambdal,par,weights){
  c<-.001
  score<-NULL
  #for (v in glob) scoreret[v]<- scoreret[v] -lambdagl*NormShiftvder(beta[v],c)
  for (v in 1:length(par)) score<- c(score,lambdal*weights[v]*NormShiftvder(par[v],c))
  return(score)
}

############################
NormShiftv <-function(x,c){
  ## used vector argument  
  y  <-sqrt(sum(x^2) + c)-sqrt(c)
  return(y)}
###############################
NormShiftvder <-function(x,c){
  ## used vector argument
  #y<-0
  y<-1 / sqrt(sum(x^2) + c)*matrix(x)
  return(y)}
####################################

Fitmixcum<-function(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda =NULL,datnew=NULL, stand=NULL){
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  if(is.null(lambda ))lambda <-c(0,0)
  
  #pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  #resp<-as.matrix(as.numeric(dat[,vresp]))
  
  ## with design if factors present
  
  dm<-GenDatSplits(dat,varloc,varmixt)
  pred<-dm$pred
  predmix<-dm$predmix
  namespred <-dm$namespred 
  namespredmix<-dm$namespredmix
  dm$locvar  ## list with placements for variables
  dm$mixvar
  listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                   mixvar=dm$namespredmix, plmixvar=dm$mixvar)
  #listlambda$locvar[[1]]
  
  resp<-as.matrix(as.numeric(dat[,vresp]))
  #####################################
  
  maxit<-500
  
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  
  
  if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
  if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
  
  ###### standardization
  if (!is.null(stand)){
    pred<-scale(pred)
    if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
  }
  ####################
  
  
  pmix<-dim(predmix)[2]
  
  #gama<-NULL
  #if(dim(predmix)[2]>1)gama<-c(0,rep(-.5,dim(predmix)[2]-1))
  
  #Xgr<-as.matrix(dat[,pred])  ### !!!!!!
  respm <- matrix(0,dim(pred)[1],k)
  for (i in 1:dim(pred)[1]) respm [i,resp[i]]<-1
  ngrgr <- rep(1,dim(pred)[1])
  
  ##start values
  
  if(sum(start)==0){ 
    gama<-0
    if(dim(predmix)[2]>1)gama<-c(0,rep(-.5,dim(predmix)[2]-1))
    thr<-(seq(0,1,1/(k-2))-.5)*.1
    parst<-c(thr,rep(0.1,p),gama) ### first thresholds
  }
  if(sum(start)!=0)parst<-start
  
  ### fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtcum, gr = NULL,respm=respm,resp=resp,
                                    pred=pred,predmix=predmix,listlambda=listlambda, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = FALSE)
  
  if (deriv ==TRUE){
    fitopt <- optim(parst, loglikmixtcum, gr = derloglikmixtcum,respm=respm,resp=resp,
                    pred=pred,predmix=predmix,listlambda=listlambda, method = "BFGS",
                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                    hessian = hessian)
    
  }
  thrtr<-fitopt$par[1:(k-1)]
  thr<-thrtr
  for(i in 2:(k-1))thr[i]<-thr[1]+sum(exp(thrtr[2:i]))
  thr
  
  #loglikmixtcum(fitopt$par,respm,resp,pred,predmix,listlambda)
  #derloglikmixtcum(fitopt$par,respm,resp,pred,predmix,listlambda)
  #par<-fitopt$par
  listn<-listlambda
  listn[[1]]<-c(0,0)
  loglikunpen<--loglikmixtcum(fitopt$par,respm,resp,pred,predmix,listn)
  loglikpen<--fitopt$value
  der<-derloglikmixtcum(fitopt$par,respm,resp,pred,predmix,listlambda)
  
  AIC<- -2*(loglikunpen-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[k:(k+p-1)]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    #hessinv<-ginv(fitopt$hessian)
    hessinv <- tryCatch(
      ginv(fitopt$hessian),
      error = function(e) {
        message("ginv failed: ", e$message)
        return(matrix(NA, nrow(fitopt$hessian), ncol(fitopt$hessian)))
      }
    )
    
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    
    #location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p],pval[1:p])
    location<-cbind( fitopt$par[k:(k+p-1)], stderr[k:(k+p-1)],zval[k:(k+p-1)],pval[k:(k+p-1)])
    #
    mixture<-cbind(fitopt$par[(p+k):length(parst)], stderr[(p+k):length(parst)],
                   zval[(p+k):length(parst)],pval[(p+k):length(parst)])
    
    
    location<-as.data.frame(location)
    
    names(location)[1] <- "Estimates"
    names(location)[2] <- "std err"
    names(location)[3] <- "z-values"
    names(location)[4] <- "p-values"
    row.names(location)<-namespred #c(seq(1,(knew-1),1),paste(pred))
    
    mixture<-as.data.frame(mixture)
    
    names(mixture)[1] <- "Estimates"
    names(mixture)[2] <- "std err"
    names(mixture)[3] <- "z-values"
    names(mixture)[4] <- "p-values"
    row.names(mixture)<-c("const",namespredmix)
  }
  
  if(hessian ==FALSE){
    #hessinv<-solve(fitopt$hessian)
    location<- fitopt$par[k:(k+p-1)]
    mixture<-fitopt$par[(p+k):length(parst)]
  }
  
  #### prediction
  loglikpred<-0
  
  if(!is.null(datnew )){
    dm<-GenDatSplits(datnew,varloc,varmixt)
    pred<-dm$pred
    predmix<-dm$predmix
    namespred <-dm$namespred 
    namespredmix<-dm$namespredmix
    dm$locvar  ## list with placements for variables
    dm$mixvar
    listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                     mixvar=dm$namespredmix, plmixvar=dm$mixvar)
    
    resp<-as.matrix(as.numeric(datnew[,vresp]))
    
    if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
    if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
    
    ###### standardization
    if (sum(lambda[[1]])>0){
      pred<-scale(pred)
      if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
    }
    ####################
    listn<-listlambda
    listn[[1]]<-c(0,0)
    respm <- matrix(0,dim(pred)[1],k)
    for (i in 1:dim(pred)[1]) respm [i,resp[i]]<-1
    ngrgr <- rep(1,dim(pred)[1])
    loglikpred<- -loglikmixtcum(fitopt$par,respm,k,pred,predmix,listn) 
    #loglikmixtcum(fitopt$par,respm,resp,pred,predmix,listn)
  }
  
  newList <- list("Loglik"= loglikunpen,"LoglikPen"= loglikpen,"AIC"=AIC, "location"=location, "mixture"= mixture,
                  "parameter"=fitopt$par, "convergence"= fitopt$convergence,
                  "der"=der,"ordered thresholds" =thr,"loglikpred"=loglikpred)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}




###################################
loglikmixtcum <- function (par,respm,resp,pred,predmix,listlambda){
  
  #### model Y <= r in predictor - x
  #### respm: responses as counts in n x k matrix 
  #### ngr:   numbers of observations at fixed value (n x 1) not needed in loglik (but score)
  
  X<-pred
  n<-dim(X)[1]
  k<- dim(respm)[2]
  
  
  ## transform thresholds
  thr<-par[1:(k-1)]
  for(i in 2:(k-1))thr[i]<-par[1]+sum(exp(par[2:i]))
  par[1:(k-1)]<-thr
  ####################
  beta<-par[1:(k-1+dim(pred)[2])]
  gama<-par[(k+dim(pred)[2]):length(par)]
  
  
  ones<- matrix(1,k-1,1)
  prob<- matrix(0,n,k)  ### cumulative
  probtot<- matrix(0,n,k)
  
  for (i in 1:n){
    etam <- cbind(diag(k-1),-ones%*%X[i,])   ### predictor - x
    prob[i,1]<- plogis(etam[1,]%*%beta)
    if(k >=3){for (j in 2:(k-1)) {
      prob[i,j]<-plogis(etam[j,]%*%beta) -plogis(etam[j-1,]%*%beta)}
    }
    prob[i,k]<- 1-sum(prob[i,])
    #ind[i,k]<-NegVal(prob[i,k])
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){probtot[i,r]<-pimix*prob[i,r]+ (1-pimix)/k}
    
    #corr<-.00000001
    #probtot[i,]<-(probtot[i,]+corr)/(1+k*corr)
  }
  #sum(probtot[i,])
  #sum(prob[i,])
  ## loglik 
  
  
  loglik <-0
  for (i in 1:n)loglik <- loglik +respm[i,]%*%log(probtot[i,])
  
  ### pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        loglik<-loglik-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*NormShiftv(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        loglik<-loglik-sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*NormShiftv(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  loglik <- -loglik ### negative log-likelihood
  #print(par)
  #print(-loglik)
  return(loglik)
}
#################################


derloglikmixtcum <- function(par, respm, resp, pred, predmix,listlambda){
  
  X <- pred
  n <- dim(X)[1]
  k <- dim(respm)[2]
  p <- dim(pred)[2]
  
  # save original par for chain rule
  par_orig <- par
  
  # transform thresholds FIRST - exactly as in log-likelihood
  thr <- par[1:(k-1)]
  for(i in 2:(k-1)) thr[i] <- par[1] + sum(exp(par[2:i]))
  par[1:(k-1)] <- thr
  
  # extract beta AFTER transformation - exactly as in log-likelihood
  beta <- par[1:(k-1+p)]
  gama <- par[(k+p):length(par)]
  
  ones <- matrix(1, k-1, 1)
  prob    <- matrix(0, n, k)
  probtot <- matrix(0, n, k)
  cumF    <- matrix(0, n, k-1)
  
  # forward pass - identical to log-likelihood
  for(i in 1:n){
    etam <- cbind(diag(k-1), -ones %*% X[i,])
    cumF[i,1] <- plogis(etam[1,] %*% beta)
    prob[i,1] <- cumF[i,1]
    if(k >= 3){
      for(j in 2:(k-1)){
        cumF[i,j] <- plogis(etam[j,] %*% beta)
        prob[i,j] <- cumF[i,j] - cumF[i,j-1]
      }
    }
    prob[i,k] <- 1 - sum(prob[i,])
    pimix <- exp(predmix[i,] %*% gama) / (1 + exp(predmix[i,] %*% gama))
    for(r in 1:k) probtot[i,r] <- pimix * prob[i,r] + (1-pimix)/k
  }
  
  # gradient accumulators w.r.t. transformed beta (length k-1+p) and gama
  derbeta <- matrix(0, k-1+p, 1)
  dergama <- matrix(0, length(gama), 1)
  
  for(i in 1:n){
    etam  <- cbind(diag(k-1), -ones %*% X[i,])
    pimix <- exp(predmix[i,] %*% gama) / (1 + exp(predmix[i,] %*% gama))
    dF    <- cumF[i,] * (1 - cumF[i,])  # length k-1
    
    for(r in 1:k){
      if(r == 1){
        dPdbeta <- dF[1] * etam[1,]
      } else if(r == k){
        dPdbeta <- -dF[k-1] * etam[k-1,]
      } else {
        dPdbeta <- dF[r] * etam[r,] - dF[r-1] * etam[r-1,]
      }
      
      weight  <- respm[i,r] * as.numeric(pimix / probtot[i,r])
      derbeta <- derbeta + weight * as.matrix(dPdbeta)
    }
    
    # gamma derivative
    for(r in 1:k){
      dergama <- dergama + respm[i,r] *
        as.numeric((prob[i,r] - 1/k) * pimix * (1-pimix) / probtot[i,r]) *
        as.matrix(predmix[i,])
    }
  }
  
  # chain rule: derbeta[1:(k-1)] is w.r.t. transformed thresholds
  # need to map back to par_orig[1:(k-1)]
  # par_orig[1]: d thr_j / d tau_1 = 1 for all j -> sum(derbeta[1:(k-1)])
  # par_orig[l] for l>=2: d thr_j / d tau_l = exp(par_orig[l]) for j>=l -> exp(par_orig[l])*sum(derbeta[l:(k-1)])
  derthr <- derbeta[1:(k-1)]
  
  derpar_orig <- matrix(0, k-1, 1)
  derpar_orig[1] <- sum(derthr)
  if(k >= 3){
    for(l in 2:(k-1)){
      derpar_orig[l] <- exp(par_orig[l]) * sum(derthr[l:(k-1)])
    }
  }
  
  # x-coefficients need no chain rule
  derbetax <- derbeta[k:(k-1+p)]
  
  der <- c(derpar_orig, derbetax, dergama)  ### - deleted
  
  ### new pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        der[listlambda$pllocvar[[v]]]<-der[listlambda$pllocvar[[v]]]-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*
          NormShiftvder(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        der[p+k-1+listlambda$plmixvar[[v]]+1]<-der[p+k-1+listlambda$plmixvar[[v]]+1]-
          sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*
          NormShiftvder(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  
  der<--der
  
  return(der)
}

######################

FitmixbinScaled <-function(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda=NULL,datnew=NULL, stand=NULL){
  #Fitmixbin <-function(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda=NULL,datnew=NULL) 
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  if(is.null(lambda ))lambda <-c(0,0)
  #pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  ### design for ordered variables
  
  dm<-GenDatSplits(dat,varloc,varmixt)
  pred<-dm$pred
  predmix<-dm$predmix
  namespred <-dm$namespred 
  namespredmix<-dm$namespredmix
  dm$locvar  ## list with placements for variables
  dm$mixvar
  listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                   mixvar=dm$namespredmix, plmixvar=dm$mixvar)
  #listlambda$locvar[[1]]
  
  
  
  resp<-as.matrix(as.numeric(dat[,vresp]))
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
  if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
  
  ###### standardization
  if (!is.null(stand)){
    pred<-scale(pred)
    if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
  }
  ####################
  maxit<-100
  
  pmix<-dim(predmix)[2]
  
  ##start values
  
  #fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = FALSE,maxit=maxit,start=0)
  #fit
  if(sum(start)==0){
    #fit<-Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit,start=0)
    gama<-NULL
    if(dim(predmix)[2]>1)gama<-c(0.5,rep(.01,dim(predmix)[2]-1))
    if(dim(predmix)[2]==1)gama<-0.5
    parst<-c(rep(0,p),1,gama,1)  ## last is scaling  ###!
  }
  
  if(sum(start)!=0)parst<-start
  
  ## fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtbinScaled, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix,
                                    listlambda=listlambda,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv ==TRUE){
    if(sum(start)==0){maxit<-100
    fitoptnd <- optim(parst, loglikmixtbinScaled, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix, 
                      listlambda=listlambda,method = "Nelder-Mead",
                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    parst<-fitoptnd$par}
    
    fitopt <- optim(parst, loglikmixtbinScaled, gr = derloglikmixtbinScaled,resp=resp,k=k,
                    pred=pred,predmix=predmix,listlambda=listlambda, method = "BFGS",
                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                    hessian = hessian)
  }
  
  #loglikmixtbinScaled(fitopt$par,resp,k,pred,predmix,listlambda)
  #derloglikmixtbinScaled(parst,resp,k,pred,predmix,listlambda)
  der<-derloglikmixtbinScaled(fitopt$par,resp,k,pred,predmix,listlambda)
  #loglikmixtbinScaled<-function(par,resp,k,pred,predmix,listlambda)
  listn<-listlambda
  listn[[1]]<-c(0,0)
  loglikunpen<--loglikmixtbinScaled(fitopt$par,resp,k,pred,predmix,listn)
  loglikpen<--fitopt$value
  
  
  AIC<- -2*(loglikunpen-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  mixture<-fitopt$par[(p+2):(length(parst)-1)]
  scaling<-fitopt$par[length(parst)]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    #hessinv <-  ginv(fitopt$hessian)
    H <- fitopt$hessian
    eps <- 1e-8
    
    hessinv<- solve(H + diag(eps, nrow(H))) 
    
    #hessinv<-ginv(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p],pval[1:p])
    mixture<-cbind(fitopt$par[(p+2):(length(parst)-1)], stderr[(p+2):(length(parst)-1)],
                   zval[(p+2):(length(parst)-1)],pval[(p+2):(length(parst)-1)])
    
    location<-as.data.frame(location)
    
    names(location)[1] <- "Estimates"
    names(location)[2] <- "std err"
    names(location)[3] <- "z-values"
    names(location)[4] <- "p-values"
    row.names(location)<-namespred 
    
    mixture<-as.data.frame(mixture)
    
    names(mixture)[1] <- "Estimates"
    names(mixture)[2] <- "std err"
    names(mixture)[3] <- "z-values"
    names(mixture)[4] <- "p-values"
    row.names(mixture)<-c("const",namespredmix)
    
    dimpar<-length(parst)
    zvalscale<-(fitopt$par[dimpar]-1)/stderr[dimpar]
    pvalscale <-(1-pnorm(abs(zvalscale), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    scaling<-cbind(fitopt$par[dimpar], stderr[dimpar],zvalscale,pvalscale)
    scaling<-as.data.frame(scaling)
    names(scaling)[1] <- "Estimates"
    names(scaling)[2] <- "std err"
    names(scaling)[3] <- "z-values"
    names(scaling)[4] <- "p-values (dev from 1)"
    
    
  }
  
  #### prediction
  loglikpred<-0
  
  if(!is.null(datnew )){
    dm<-GenDatSplits(datnew,varloc,varmixt)
    pred<-dm$pred
    predmix<-dm$predmix
    namespred <-dm$namespred 
    namespredmix<-dm$namespredmix
    dm$locvar  ## list with placements for variables
    dm$mixvar
    listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                     mixvar=dm$namespredmix, plmixvar=dm$mixvar)
    
    resp<-as.matrix(as.numeric(datnew[,vresp]))
    
    if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
    if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
    
    ###### standardization
    if (sum(lambda[[1]])>0){
      pred<-scale(pred)
      if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
    }
    ####################
    listn<-listlambda
    listn[[1]]<-c(0,0)
    
    loglikpred<- -loglikmixtbinScaled(fitopt$par,resp,k,pred,predmix,listn) 
    
  }
  
  
  newList <- list("Loglik"= loglikunpen,"LoglikPen"= loglikpen,"AIC"=AIC, "location"=location, "mixture"= mixture,"scaling"=scaling,
                  "parameter"=fitopt$par, "convergence"= fitopt$convergence,"der"=der,"loglikpred"=loglikpred)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

#################################

loglikmixtbinScaledUnrestricted<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # scaling parameter not restricted
  
  
  #zgamma<- 1  # not yet needed'  
  ## par (beta, intercept,gamma(mixture),scale)
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  scalepar<-par[length(par)]
  
  ### thresholds
  thr<-NULL
  for (r in 2:k) thr<-c(thr,scalepar*log((k-r+1)/(r-1)))
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds  with 0
  gama<-par[(p+2):(length(par)-1)]   ### -1
  
  #dim(predmix)
  #length(par)
  ###### probabilities
  
  n<-dim(pred)[1]
  
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){prob[i,r]<-pimix*expterm[1,r]/sumn+ (1-pimix)/k}
    #sum(prob[i,])
    
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ### pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        loglik<-loglik-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*NormShiftv(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        loglik<-loglik-sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*NormShiftv(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  #newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(loglikneg)}
#####################################################

derloglikmixtbinScaledUnrestricted<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # scaling parameter not restricted
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  scalepar<-par[length(par)]
  
  ### thresholds
  thr<-NULL
  thrunscaled<-0  ## first is 0
  for (r in 2:k) {thr<-c(thr,scalepar*log((k-r+1)/(r-1)))
  thrunscaled<-c(thrunscaled,log((k-r+1)/(r-1)))}
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds with 0  
  gama<-par[(p+2):(length(par)-1)]
  
  
  n<-dim(pred)[1]
  
  prob<- matrix(0,n,k)
  mixtprob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  #der2<-matrix(0,1,1)
  der1mixt<-matrix(0,p+1,1)  # for beta, beta_0
  der2mixt<-matrix(0,dim(predmix)[2],1) # for gamma
  der3mixt<-0
  
  der<-matrix(0,p+1+dim(predmix)[2],1)
  #dim #
  #loglik<-0
  
  for(i in 1:n){
    
    ###### probabilities adj categories
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sumn
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    sum(prob[i,])
    
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der1 <- as.matrix(((cat-1) - sum1) * pred[i,])
    
    ### thresholds
    
    #for (j in 2:k){
    #  ind1<-0
    #  if(cat>=j)ind1<-1
    #  der2[j-1]<-ind1-sum(prob[i,j:k])    }
    
    #sumdum<-0
    #for (j in 2:k)sumdum<-sumdum+prob[i,j]*(j-1)
    #der2<-cat-1-sumdum
    
    der2 <- 0
    for (j in 2:k){
      ind1 <- 0
      if(cat >= j) ind1 <- 1
      der2 <- der2 + ind1 - sum(prob[i, j:k])
    }
    
    
    deradj<-c(der1,der2)
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){mixtprob[i,r]<-pimix*prob[i,r]+ (1-pimix)/k}
    mixtsel<-mixtprob[i,][resp[i]]
    #sum(prob[i,])
    #sum(mixtprob[i,])
    der1mixt<-der1mixt+as.numeric(pimix*prob[i,][resp[i]])*matrix(deradj)/mixtsel
    
    der2mixt<-der2mixt+as.numeric((prob[i,][resp[i]]-1/k)*pimix*(1-pimix))*matrix(predmix[i,])/
      (mixtsel)
    
    ###scaling
    
    der3mixtloc<-sum(thrunscaled[1:(resp[i])])
    for (s in 2:k)der3mixtloc<-der3mixtloc-(prob[i,][s])*sum(thrunscaled[1:s])
    der3mixtloc<-  pimix*(prob[i,][resp[i]])*der3mixtloc/mixtsel
    der3mixt<-der3mixt+der3mixtloc
    
  }## end i
  
  #####  
  der<-c(der1mixt,der2mixt,der3mixt)
  
  ### new pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        der[listlambda$pllocvar[[v]]]<-der[listlambda$pllocvar[[v]]]-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*
          NormShiftvder(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        der[p+1+listlambda$plmixvar[[v]]+1]<-der[p+1+listlambda$plmixvar[[v]]+1]-
          sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*
          NormShiftvder(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  
  
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}




loglikmixtbinScaled<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  ## par (beta, intercept,gamma(mixture),scale)
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  scalepar<-exp(par[length(par)])   ### now exp!
  
  ### thresholds
  thr<-NULL
  for (r in 2:k) thr<-c(thr,scalepar*log((k-r+1)/(r-1)))
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds  with 0
  gama<-par[(p+2):(length(par)-1)]   ### -1
  
  #dim(predmix)
  #length(par)
  ###### probabilities
  
  n<-dim(pred)[1]
  
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){prob[i,r]<-pimix*expterm[1,r]/sumn+ (1-pimix)/k}
    #sum(prob[i,])
    
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  ### pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        loglik<-loglik-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*NormShiftv(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        loglik<-loglik-sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*NormShiftv(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  #newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(loglikneg)}
#####################################################

derloglikmixtbinScaled<-function(par,resp,k,pred,predmix,listlambda){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  scalepar<-exp(par[length(par)])   ### now exp!
  
  ### thresholds
  thr<-NULL
  thrunscaled<-0  ## first is 0
  for (r in 2:k) {thr<-c(thr,scalepar*log((k-r+1)/(r-1)))
  thrunscaled<-c(thrunscaled,log((k-r+1)/(r-1)))}
  thr<-thr+par[p+1]
  beta0<-c(0,thr)   ### thresholds with 0  
  gama<-par[(p+2):(length(par)-1)]
  
  
  n<-dim(pred)[1]
  
  prob<- matrix(0,n,k)
  mixtprob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  #der2<-matrix(0,1,1)
  der1mixt<-matrix(0,p+1,1)  # for beta, beta_0
  der2mixt<-matrix(0,dim(predmix)[2],1) # for gamma
  der3mixt<-0
  
  der<-matrix(0,p+1+dim(predmix)[2],1)
  #dim #
  #loglik<-0
  
  for(i in 1:n){
    
    ###### probabilities adj categories
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
    sumn <-1
    for (r in 2:k){sumn<-sumn+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sumn
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    sum(prob[i,])
    
    #loglik<-loglik+log(prob[i,])[resp[i]]
    cat<-resp[i]
    
    ### covariates
    der1<-(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der1 <- as.matrix(((cat-1) - sum1) * pred[i,])
    
    ### thresholds
    
    #for (j in 2:k){
    #  ind1<-0
    #  if(cat>=j)ind1<-1
    #  der2[j-1]<-ind1-sum(prob[i,j:k])    }
    
    #sumdum<-0
    #for (j in 2:k)sumdum<-sumdum+prob[i,j]*(j-1)
    #der2<-cat-1-sumdum
    
    der2 <- 0
    for (j in 2:k){
      ind1 <- 0
      if(cat >= j) ind1 <- 1
      der2 <- der2 + ind1 - sum(prob[i, j:k])
    }
    
    
    deradj<-c(der1,der2)
    
    pimix<-exp(predmix[i,]%*%gama)/(1+exp(predmix[i,]%*%gama))
    
    for (r in 1:k){mixtprob[i,r]<-pimix*prob[i,r]+ (1-pimix)/k}
    mixtsel<-mixtprob[i,][resp[i]]
    #sum(prob[i,])
    #sum(mixtprob[i,])
    der1mixt<-der1mixt+as.numeric(pimix*prob[i,][resp[i]])*matrix(deradj)/mixtsel
    
    der2mixt<-der2mixt+as.numeric((prob[i,][resp[i]]-1/k)*pimix*(1-pimix))*matrix(predmix[i,])/
      (mixtsel)
    
    ###scaling
    
    der3mixtloc<-sum(thrunscaled[1:(resp[i])])*scalepar
    for (s in 2:k)der3mixtloc<-der3mixtloc-(prob[i,][s])*sum(thrunscaled[1:s])*scalepar
    der3mixtloc<-  pimix*(prob[i,][resp[i]])*der3mixtloc/mixtsel
    der3mixt<-der3mixt+der3mixtloc
    #der3mixt<-der3mixt*scalepar  ###!
  }## end i
  
  #####  
  der<-c(der1mixt,der2mixt,der3mixt)
  
  ### new pen  ## without weights
  c<-.001
  if (listlambda$lam[[1]] >0){
    lambda<-listlambda$lam
    varn<-length(listlambda$pllocvar)   
    if (varn >0){
      for (v in 1:varn)  
        der[listlambda$pllocvar[[v]]]<-der[listlambda$pllocvar[[v]]]-sqrt(length(listlambda$pllocvar[[v]]))*lambda[1]*
          NormShiftvder(beta[listlambda$pllocvar[[v]]],c)
    }
  }
  if (listlambda$lam[[2]] >0){     
    lambda<-listlambda$lam
    if ((varn >0)&(length(gama)>1)){
      for (v in 1:length(listlambda$plmixvar))  
        der[p+1+listlambda$plmixvar[[v]]+1]<-der[p+1+listlambda$plmixvar[[v]]+1]-
          sqrt(length(listlambda$plmixvar[[v]]))*lambda[2]*
          NormShiftvder(gama[listlambda$plmixvar[[v]]+1],c)
    }
  }
  
  
  
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}



####check
















#########################
#You can verify with:
#library(numDeriv)

#gr_num  <- grad(loglikmixtcum, par_start, respm=respm, resp=resp, pred=pred, predmix=predmix)
#gr_anal <- derloglikmixtcum(par_start, respm, resp, pred, predmix)
#cbind(gr_num, gr_anal)

########################### old
FitmixbinOld <-function(vresp,varloc,varmixt,dat,k,deriv,hessian,start,lambda=NULL,datnew=NULL){
  
  ### vresp:responses variable
  ### varloc: predictors location 
  ### varmixt: predictors mixture 
  ### k: number of categories, 1,2,...,k
  ### dat: data
  ### deriv: if 'der' derivatives used
  ### hessian: TRUE or FALSE
  ### start: start values, none if 0
  
  disp<-0  ## obsolete value
  if(is.null(lambda ))lambda <-c(0,0)
  #pred <-as.matrix(dat[,varloc])
  #predmix<-as.matrix(dat[,varmixt])
  #if ( is.null(varmixt))predmix<-NULL
  ### design for ordered variables
  
  dm<-GenDatSplits(dat,varloc,varmixt)
  pred<-dm$pred
  predmix<-dm$predmix
  namespred <-dm$namespred 
  namespredmix<-dm$namespredmix
  dm$locvar  ## list with placements for variables
  dm$mixvar
  listlambda<-list(lam=lambda,locvar=dm$namespred,pllocvar=dm$locvar,
                   mixvar=dm$namespredmix, plmixvar=dm$mixvar)
  #listlambda$locvar[[1]]
  
  
  
  resp<-as.matrix(as.numeric(dat[,vresp]))
  p<-dim(pred)[2]
  #beta<-rep(.01,p)  ### first are covariate weights
  #beta0<-rep(.02,(k-1)) 
  if ( !is.null(predmix))predmix<-as.matrix(cbind(rep(1,dim(pred)[1]),predmix))
  if ( is.null(predmix))predmix<-as.matrix(rep(1,dim(pred)[1]))
  
  ###### standardization
  if (sum(lambda[[1]])>0){
    pred<-scale(pred)
    if ( !is.null(varmixt))predmix[,2:dim(predmix)[2]]<-scale(predmix[,2:dim(predmix)[2]])
  }
  ####################
  maxit<-500
  
  pmix<-dim(predmix)[2]
  
  ##start values
  
  #fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = FALSE,maxit=maxit,start=0)
  #fit
  if(sum(start)==0){
    fit<-Fitadj(vresp,varloc,dat,k,deriv,hessian,maxit,start=0)
    fit
    
    gama<-NULL
    if(dim(predmix)[2]>1)gama<-c(0,rep(.01,dim(predmix)[2]-1))
    if(dim(predmix)[2]==1)gama<-0.1
    parst<-c(fit$parameter[1:p]/2,0,gama)
  }
  
  if(sum(start)!=0)parst<-start
  
  ## fit
  
  if (deriv ==FALSE)fitopt <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix,
                                    listlambda=listlambda,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv ==TRUE){
    if(sum(start)==0){maxit<-100
    fitoptnd <- optim(parst, loglikmixtbin, gr = NULL,resp=resp,k=k,pred=pred,predmix=predmix, 
                      listlambda=listlambda,method = "Nelder-Mead",
                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    parst<-fitoptnd$par}
    
    fitopt <- optim(parst, loglikmixtbin, gr = derloglikmixtbin,resp=resp,k=k,
                    pred=pred,predmix=predmix,listlambda=listlambda, method = "BFGS",
                    lower = -Inf, upper = Inf,control = list(maxit=maxit), 
                    hessian = hessian)
  }
  
  #loglikmixtadjcat(parst,resp,k,pred,predmix)
  der<-derloglikmixtbin(fitopt$par,resp,k,pred,predmix,listlambda)
  #derloglikmixtbin(parst,resp,k,pred,predmix)
  listn<-listlambda
  listn[[1]]<-c(0,0)
  loglikunpen<--loglikmixtbin(fitopt$par,resp,k,pred,predmix,listn)
  loglikpen<--fitopt$value
  
  
  AIC<- -2*(loglikunpen-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  mixture<-fitopt$par[(p+2):length(parst)]
  
  if(hessian ==TRUE){
    #hessinv<-solve(fitopt$hessian)
    #hessinv <-  ginv(fitopt$hessian)
    H <- fitopt$hessian
    eps <- 1e-8
    
    hessinv<- solve(H + diag(eps, nrow(H))) 
    
    #hessinv<-ginv(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std 
    zval<-fitopt$par/stderr
    pval <-(1-pnorm(abs(zval), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE))*2
    
    location<-cbind( fitopt$par[1:p], stderr[1:p],zval[1:p],pval[1:p])
    mixture<-cbind(fitopt$par[(p+2):length(parst)], stderr[(p+2):length(parst)],
                   zval[(p+2):length(parst)],pval[(p+2):length(parst)])
    
    location<-as.data.frame(location)
    
    names(location)[1] <- "Estimates"
    names(location)[2] <- "std err"
    names(location)[3] <- "z-values"
    names(location)[4] <- "p-values"
    row.names(location)<-varloc 
    
    mixture<-as.data.frame(mixture)
    
    names(mixture)[1] <- "Estimates"
    names(mixture)[2] <- "std err"
    names(mixture)[3] <- "z-values"
    names(mixture)[4] <- "p-values"
    row.names(mixture)<-c("const",namespredmix)
    
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= loglikunpen,"LoglikPen"= loglikpen,"AIC"=AIC, "location"=location, "mixture"= mixture,
                  "parameter"=fitopt$par, "convergence"= fitopt$convergence,"der"=der)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

#################################


