###############################################################################################################################
####    THIS PROGRAM PROVIDES FUNCTIONS TO MAKE PREDICTION UNDER THE SKEW-NORMAL SAE MODEL                                #####
#####                                                                                                                     #####
#####                                                                                                                     #####
#####   Autheur: Mamadou S. Diallo                                                                                        #####
#####   Date:    October 2013                                                                                             #####
###############################################################################################################################



### Loading the libraries needed for this program ###

library(Matrix)
library(MASS)
library(mvtnorm)
library(tmvtnorm)
#library(tmg)
library(numDeriv)




## This function computes the best predictor when errors are SN 
CSN_Predict_d <- function(d, L, samp, beeta=beeta, sig_e=sig_e, lambda_e=lambda_e, gamma_nd=gamma_nd, mu_u=mu_u, mu_e=mu_e, 
                          alpha=alpha, beta=beta, B0=B0, nu=nu, gamma_inv=gamma_inv, PovLine=12){

  pop    <- samp[samp$Area%in%d,]
  pop.yr <- as.matrix(pop[pop$Sel%in%0,])
  pop.ys <- as.matrix(pop[pop$Sel%in%1,])
  y_ds   <- matrix(pop.ys[,1])
  xbetar <- (pop.yr[,3:5]%*%beeta)[,1]

  Ndr    <- nrow(pop.yr)
  Nds    <- nrow(pop.ys)
  res_d  <- y_ds - ((pop.ys[,3:5]%*%beeta)[,1] + mu_u + mu_e)
  mu_drs <- xbetar + (mu_u + mu_e + gamma_nd*sum(res_d))*rep(1,Ndr)
  nu_0s  <- -as.vector(nu%*%res_d)
  
  V1     <- replicate(L,rnorm(n=Ndr,mean=0,sd=sqrt(alpha))+rnorm(n=1,mean=0,sd=sqrt(beta)))
  V0r    <- abs(replicate(L,rnorm(n=Ndr,mean=0,sd=sqrt(1+lambda_e^2))))
  #V0s    <- t(rtmg(n=L, M=gamma_inv, r=rep(0,Nds+1), initial=nu_0s+0.1, f = diag(Nds+1), g = -nu_0s, burn.in = 50))
  V0s    <- t(rtmvnorm(n=L, mean=rep(0,param$nd+1), H=gamma_inv, lower=nu_0s, algorithm="gibbs", burn.in.samples=500, thinning=5))
  V0     <- rbind(V0r,V0s)

  y_dr   <- mu_drs + B0%*%V0  + V1

  y_pred <- rbind(y_dr, matrix(y_ds,Nds,L))

  F_d     <- matrix(data=0, nrow=1, ncol=3)

  ### Predicted values of Y
  F_d[1,1] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^0)))  #poverty incidence
  F_d[1,2] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^1)))  #poverty gap
  F_d[1,3] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^2)))  #poverty severity

  ### Data to return
  return(F_d)
}

CSN_Predict <- function(samp, beeta, sig_e, sig_u, lambda_u, lambda_e, adjust=TRUE){

  #Computes input values for CSN_Predict
  delta_u  <- lambda_u/sqrt(1+lambda_u^2)
  delta_e  <- lambda_e/sqrt(1+lambda_e^2)
  gamma_nd <- sig_u^2/(sig_e^2+param$nd*sig_u^2)
  alpha    <- sig_e^2*(1-delta_e^2)
  den      <- (1+param$nd*gamma_nd*lambda_e^2)*sig_u^2+gamma_nd*(lambda_u*sig_e)^2
  beta     <- gamma_nd*(sig_u*sig_e)^2/den
  b01      <- -gamma_nd*lambda_e*sig_e*sig_u^2/den
  b02      <-  gamma_nd*lambda_u*sig_u*sig_e^2/den
  B0       <- cbind((lambda_e*sig_e/(1+lambda_e^2))*diag(param$Nd-param$nd),b01*matrix(1,param$Nd-param$nd,param$nd),b02*rep(1,param$Nd-param$nd))
  gamma_0s <- matrix(rbind(cbind(diag(param$nd)+lambda_e^2*gamma_nd,-lambda_e*lambda_u*(sig_e/sig_u)*gamma_nd*rep(1,param$nd)),
              c(-lambda_e*lambda_u*(sig_e/sig_u)*gamma_nd*rep(1,param$nd),1+lambda_u^2*(sig_e/sig_u)^2*gamma_nd)),param$nd+1,param$nd+1)
  gamma_inv<- matrix(rbind(cbind(diag(param$nd)-(gamma_nd*lambda_e^2*sig_u^2/den)*matrix(1,param$nd,param$nd),
                     (lambda_e*lambda_u*sig_e*sig_u*gamma_nd/den)*rep(1,param$nd)),c((lambda_e*lambda_u*sig_e*sig_u*gamma_nd/den)*rep(1,param$nd),
                     (1+param$nd*gamma_nd*lambda_e^2)*sig_u^2/den)),param$nd+1,param$nd+1)

  #Adjust the mean to get zero-mean error model
  if (adjust == "TRUE") {
    mu_u <- -delta_u*sig_u*sqrt(2/pi)
    mu_e <- -delta_e*sig_e*sqrt(2/pi)
  } else {
    mu_u <- 0
    mu_e <- 0
  }
  nu      <- rbind((lambda_e/sig_e)*(diag(param$nd)-gamma_nd*matrix(1,param$nd,param$nd)),(lambda_u/sig_u)*gamma_nd)

  #Computes the predictor for all the small  areas
  y_pred <- t(sapply(1:param$d,FUN=CSN_Predict_d, L=param$L, beeta=beeta, samp=samp, sig_e=sig_e, lambda_e=lambda_e, gamma_nd=gamma_nd, mu_u=mu_u, mu_e=mu_e, 
                          alpha=alpha, beta=beta, B0=B0, nu=nu, gamma_inv=gamma_inv, PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(y_pred)
  names(obj_pred) <- c("I","G","S")
  return(obj_pred)
  
}


## This function computes the best predictor using the conditional approach when errors are SN 

### The function Phi.d computes an approximation of Phi for small are d 
### Reference: volume1: Models and Applications 2nd ed. author: Samuel Kotz, N. Balakrishnan and N. L. Jonhson
Phi.d <- function(u0,zj,cj){
  phi.j <- function(u0,zj,cj) dnorm(u0)*prod(pnorm((zj-cj*u0)/sqrt(1-cj^2)))
  return(mapply(phi.j, u0, MoreArgs=list(zj=zj,cj=cj), SIMPLIFY=TRUE))
}

CSN_Cond_Predict_d <- function(d,L,samp,beeta,sig_u,lambda_u,sig_e,lambda_e,nu,gamma1,gamma2,Gamma,cj,adjust=TRUE,PovLine=12){

  pop    <- samp[samp$Area%in%d,]
  pop.yr <- as.matrix(pop[pop$Sel%in%0,])
  pop.ys <- as.matrix(pop[pop$Sel%in%1,])
  y_ds   <- as.vector(pop.ys[,1])
  x_ds   <- as.vector((pop.ys[,3:5]%*%beeta)[,1])
  xbetar <- as.vector((pop.yr[,3:5]%*%beeta)[,1])

  Ndr     <- nrow(pop.yr)
  Nds     <- nrow(pop.ys)
       
  #Adjust the mean to get zero-mean error model
  delta_u <- lambda_u/sqrt(1+lambda_u^2)
  delta_e <- lambda_e/sqrt(1+lambda_e^2)
  if (adjust == "TRUE") {
    mu_u <- -delta_u*sig_u*sqrt(2/pi)
    mu_e <- -delta_e*sig_e*sqrt(2/pi)
  } else {
    mu_u <- 0
    mu_e <- 0
  }
  
  nu       <- nu%*%(y_ds-(x_ds+mu_u+mu_e))
  
  Phi.dt <- function(t,cj,nu,gamma1,Gamma){
    V  <- diag(diag(Gamma))
    Da <- c(lambda_u/sig_u,-rep(lambda_e/sig_e,param$nd))
    zj <- as.vector(t(Da*sig_e^2*gamma1*t-nu)%*%diag(diag(V)^(-1/2)))
    #limit and nd.points were chosen to give an excellent approximation more than 20 digits after the dot
    limit      <- 25
    nb.points  <- 1000
    int.length <- 2*limit/nb.points
    start.int  <- -limit + int.length*0.5
    end.int    <-  limit - int.length*0.5
    univers    <- seq(from=start.int,to=end.int,by=int.length)
    return(int.length*sum(Phi.d(u0=univers,zj=zj,cj=cj)))
  }

  ## Area effects
  res_BLP <- as.vector(y_ds-x_ds)
  res_BP  <- as.vector(y_ds-(x_ds+mu_u+mu_e))
  u_EBLP  <- as.numeric(gamma2*sum(res_BLP))  #as.numeric(gamma2*sum(res_BLP[pop.ys[,2]==d])) 

  dev.Phi.d <- grad(func=Phi.dt, x=0,cj=cj,nu=nu,gamma1=gamma1,Gamma=Gamma)
  est.Phi.d <- Phi.dt(t=0,cj=cj,nu=nu,gamma1=gamma1,Gamma=Gamma)
  if (est.Phi.d != 0) {
    u_EBP <- sum(mu_u, as.numeric(gamma1*sum(res_BP)),dev.Phi.d/est.Phi.d, na.rm=TRUE) #sum(mu_u, as.numeric(gamma1*sum(res_BP[pop.ys[,2]==d])),dev.Phi.d/est.Phi.d, na.rm=TRUE)
  } else {
    u_EBP <- mu_u + as.numeric(gamma1*sum(res_BP)) #mu_u + as.numeric(gamma1*sum(res_BP[pop.ys[,2]==d]))
  }

  ## Unit level errors 
  e_dr   <- sig_e*replicate(L,delta_e*abs(rnorm(n=Ndr))+sqrt(1-delta_e^2)*rnorm(n=Ndr))

  mu_dr     <- as.vector(xbetar + mu_e)
  y_dr_EBLP <- mu_dr + u_EBLP + e_dr
  y_dr_EBP  <- mu_dr + u_EBP  + e_dr

  y_EBLP <- rbind(y_dr_EBLP, matrix(y_ds,Nds,L))
  y_EBP  <- rbind(y_dr_EBP, matrix(y_ds,Nds,L))

  F_d     <- matrix(data=0, nrow=1, ncol=6)

  ### Predicted values of Y
  F_d[1,1] <- mean(colMeans((exp(y_EBP)<PovLine)*(((PovLine-exp(y_EBP))/PovLine)^0)))  #poverty incidence
  F_d[1,2] <- mean(colMeans((exp(y_EBP)<PovLine)*(((PovLine-exp(y_EBP))/PovLine)^1)))  #poverty gap
  F_d[1,3] <- mean(colMeans((exp(y_EBP)<PovLine)*(((PovLine-exp(y_EBP))/PovLine)^2)))  #poverty severity
  F_d[1,4] <- mean(colMeans((exp(y_EBLP)<PovLine)*(((PovLine-exp(y_EBLP))/PovLine)^0)))  #poverty incidence
  F_d[1,5] <- mean(colMeans((exp(y_EBLP)<PovLine)*(((PovLine-exp(y_EBLP))/PovLine)^1)))  #poverty gap
  F_d[1,6] <- mean(colMeans((exp(y_EBLP)<PovLine)*(((PovLine-exp(y_EBLP))/PovLine)^2)))  #poverty severity

  ### Data to return
  return(F_d)
}


CSN_Cond_Predict <- function(samp, beeta, sig_u, lambda_u, sig_e, lambda_e, adjust=TRUE, PovLine=12){

  I      <- diag(param$nd)
  J      <- matrix(1,param$nd,param$nd)
  delta_u  <- lambda_u/sqrt(1+lambda_u^2)
  delta_e  <- lambda_e/sqrt(1+lambda_e^2)
  gamma1 <- sig_u^2/(sig_e^2+param$nd*sig_u^2)
  gamma2 <- sig_u^2*(1-2*delta_u^2/pi)/(sig_e^2*(1-2*delta_e^2/pi)+param$nd*sig_u^2*(1-2*delta_u^2/pi))
  nu     <- -rbind((lambda_u/sig_u)*gamma1*t(rep(1,param$nd)),(lambda_e/sig_e)*(I-gamma1*J))
  Gamma  <- matrix(cbind(c(1+lambda_u^2*(1-param$nd*gamma1),-lambda_u*lambda_e*(sig_e/sig_u)*gamma1*rep(1,param$nd)),
                 rbind(-lambda_u*lambda_e*(sig_e/sig_u)*gamma1*t(rep(1,param$nd)),I+lambda_e^2*gamma1*J)),param$nd+1,param$nd+1)
  R       <- diag(1/sqrt(diag(Gamma)))%*%Gamma%*%diag(1/sqrt(diag(Gamma)))
  cj      <- c(-lambda_u*sig_e/sqrt(sig_e^2+param$nd*sig_u^2+(lambda_u*sig_e)^2),rep(lambda_e*sig_u/sqrt(sig_e^2+param$nd*sig_u^2+(lambda_e*sig_u)^2),param$nd))
  
  F_pred  <- t(sapply(1:param$d,CSN_Cond_Predict_d,L=param$L,samp=samp,beeta=beeta,sig_u=sig_u,lambda_u=lambda_u,sig_e=sig_e,lambda_e=lambda_e,
                       nu=nu,gamma1=gamma1,gamma2=gamma2,Gamma=Gamma,cj=cj,adjust=TRUE,PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(F_pred)
  names(obj_pred) <- c("I_EBP","G_EBP","S_EBP","I_EBLP","G_EBLP","S_EBLP")
  return(obj_pred)

}


## This function computes the best predictor when errors are SN 
SN_Predict_d <- function(d, L, samp, beeta=beeta, sig_e=sig_e, gamma_nd=gamma_nd, mu_u=mu_u, mu_e=mu_e, alpha=alpha, beta=beta, b_u=b_u, b_e=b_e, PovLine=12){

  pop     <- samp[samp$Area%in%d,]
  pop.yr  <- as.matrix(pop[pop$Sel%in%0,])
  pop.ys  <- as.matrix(pop[pop$Sel%in%1,])
  y_ds    <- as.vector(matrix(pop.ys[,1]))
  xbetar  <- as.vector((pop.yr[,3:5]%*%beeta)[,1])
  Ndr     <- nrow(pop.yr)
  Nds     <- nrow(pop.ys)

  SN_Predict_d_single <- function(Ndr, Nds, pop.ys, y_ds, xbetar, L, beeta=beeta, sig_e=sig_e, gamma_nd=gamma_nd, mu_u=mu_u, mu_e=mu_e,
                             alpha=alpha, beta=beta, b_u=b_u, b_e=b_e, PovLine=PovLine){

    res_d   <- y_ds - (as.vector((pop.ys[,3:5]%*%beeta)[,1]) + mu_u + mu_e + b_u*matrix(replicate(L,abs(rnorm(n=1,mean=0,sd=1))),Nds,L,byrow=TRUE) + 
                       replicate(L, b_e*abs(rnorm(n=Nds,mean=0,sd=1))))
    mu_drs  <- xbetar + mu_u + mu_e + b_u*matrix(replicate(L,abs(rnorm(n=1,mean=0,sd=1))),Ndr,L,byrow=TRUE) + b_e*replicate(L,abs(rnorm(n=Ndr,mean=0,sd=1))) + 
               matrix(gamma_nd*colSums(res_d),Ndr,L)
    V1      <- replicate(L, rnorm(n=Ndr,mean=0,sd=sqrt(alpha))+rnorm(n=1,mean=0,sd=sqrt(beta)))
    y_dr    <- mu_drs + V1
    y_pred  <- rbind(y_dr, matrix(y_ds,Nds,L))
    F_d     <- matrix(data=0, nrow=1, ncol=3)
    F_d[1,1] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^2)))  #poverty severity
    return(as.vector(F_d))
  
  }

  f_pred_single <- SN_Predict_d_single(Ndr=Ndr,Nds=Nds,pop.ys=pop.ys,y_ds=y_ds,xbetar=xbetar,L=L,beeta=beeta,sig_e=sig_e,gamma_nd=gamma_nd, 
                                                mu_u=mu_u,mu_e=mu_e,alpha=alpha,beta=beta,b_u=b_u,b_e=b_e,PovLine=PovLine)

  SN_Predict_d_double <- function(Ndr, Nds, pop.ys, y_ds, xbetar, L, beeta=beeta, sig_e=sig_e, gamma_nd=gamma_nd, mu_u=mu_u, mu_e=mu_e,
                             alpha=alpha, beta=beta, b_u=b_u, b_e=b_e, PovLine=PovLine){

    res_d   <- y_ds - (as.vector((pop.ys[,3:5]%*%beeta)[,1]) + mu_u + mu_e + b_u*abs(rnorm(n=1,mean=0,sd=1)) + b_e*abs(rnorm(n=Nds,mean=0,sd=1)))
    mu_drs  <- xbetar + mu_u + mu_e + b_u*abs(rnorm(n=1,mean=0,sd=1)) + b_e*abs(rnorm(n=Ndr,mean=0,sd=1)) + rep(gamma_nd*sum(res_d),Ndr)
    V1      <- replicate(L, rnorm(n=Ndr,mean=0,sd=sqrt(alpha))+rnorm(n=1,mean=0,sd=sqrt(beta)))
    y_dr    <- mu_drs + V1
    y_pred  <- rbind(y_dr, matrix(y_ds,Nds,L))
    F_d     <- matrix(data=0, nrow=1, ncol=3)
    F_d[1,1] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^2)))  #poverty severity
    return(as.vector(F_d))
  
  }

  f_pred_double <- rowMeans(replicate(L,SN_Predict_d_double(Ndr=Ndr,Nds=Nds,pop.ys=pop.ys,y_ds=y_ds,xbetar=xbetar,L=L,beeta=beeta,sig_e=sig_e,gamma_nd=gamma_nd, 
                                                mu_u=mu_u,mu_e=mu_e,alpha=alpha,beta=beta,b_u=b_u,b_e=b_e,PovLine=12)))

  ### Predicted values of Y
  F_d      <- matrix(data=0, nrow=1, ncol=6)
  F_d[1,1] <- f_pred_single[1]                                                           #poverty incidence
  F_d[1,2] <- f_pred_single[2]                                                           #poverty gap
  F_d[1,3] <- f_pred_single[3]                                                           #poverty severity
  F_d[1,4] <- f_pred_double[1]                                                           #poverty incidence
  F_d[1,5] <- f_pred_double[2]                                                           #poverty gap
  F_d[1,6] <- f_pred_double[3]                                                           #poverty severity

  ### Data to return
  return(F_d)
}


SN_Predict <- function(samp, beeta, sig_e, sig_u, lambda_u, lambda_e, adjust=TRUE){

  #Computes input values for NM_Predict
  delta_u  <- lambda_u/sqrt(1+lambda_u^2)
  delta_e  <- lambda_e/sqrt(1+lambda_e^2)
  gamma_nd <- sig_u^2*(1-delta_u^2)/(sig_e^2*(1-delta_e^2)+param$nd*sig_u^2*(1-delta_u^2))
  alpha    <- sig_e^2*(1-delta_e^2)
  beta     <- alpha*gamma_nd
  b_u      <- delta_u*sig_u
  b_e      <- delta_e*sig_e

  #Adjust the mean to get zero-mean error model
  if (adjust == "TRUE") {
    mu_u <- -delta_u*sig_u*sqrt(2/pi)
    mu_e <- -delta_e*sig_e*sqrt(2/pi)
  } else {
    mu_u <- 0
    mu_e <- 0
  }

  #Computes the predictor for all the small  areas
  y_pred <- t(sapply(1:param$d,FUN=SN_Predict_d,L=param$L,samp=samp,beeta=beeta,sig_e=sig_e,gamma_nd=gamma_nd,mu_u=mu_u,mu_e=mu_e,
                               alpha=alpha,beta=beta,b_u=b_u, b_e=b_e,PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(y_pred)
  names(obj_pred) <- c("I1","G1","S1","I2","G2","S2")
  return(obj_pred)
  
}



## This function computes the best predictor when errors are normal
NM_Predict_d <- function(d, L, samp, beeta=beeta, sig_e=sig_e, gamma_nd=gamma_nd, alpha=alpha, beta=beta, PovLine=12){

  pop     <- samp[samp$Area%in%d,]
  pop.yr  <- as.matrix(pop[pop$Sel%in%0,])
  pop.ys  <- as.matrix(pop[pop$Sel%in%1,])
  y_ds    <- matrix(pop.ys[,1])
  xbetar  <- (pop.yr[,3:5]%*%beeta)[,1]

  Ndr     <- nrow(pop.yr)
  Nds     <- nrow(pop.ys)
  res_d   <- y_ds - (pop.ys[,3:5]%*%beeta)[,1]
  mu_drs  <- xbetar + gamma_nd*sum(res_d)*rep(1,Ndr)
  
  V1      <- replicate(L,rnorm(n=Ndr,mean=0,sd=sqrt(alpha)) + rnorm(n=1,mean=0,sd=sqrt(beta)))

  y_dr    <- mu_drs + V1
  y_pred  <- rbind(y_dr, matrix(y_ds,Nds,L))

  F_d     <- matrix(data=0, nrow=1, ncol=9)

  ### Predicted values of Y
  F_d[1,1] <- mean((exp(pop[,1])<PovLine)*(((PovLine-exp(pop[,1]))/PovLine)^0))          #poverty incidence
  F_d[1,2] <- mean((exp(pop[,1])<PovLine)*(((PovLine-exp(pop[,1]))/PovLine)^1))          #poverty gap
  F_d[1,3] <- mean((exp(pop[,1])<PovLine)*(((PovLine-exp(pop[,1]))/PovLine)^2))          #poverty severity
  F_d[1,4] <- mean((exp(pop.ys[,1])<PovLine)*(((PovLine-exp(pop.ys[,1]))/PovLine)^0))    #poverty incidence
  F_d[1,5] <- mean((exp(pop.ys[,1])<PovLine)*(((PovLine-exp(pop.ys[,1]))/PovLine)^1))    #poverty gap
  F_d[1,6] <- mean((exp(pop.ys[,1])<PovLine)*(((PovLine-exp(pop.ys[,1]))/PovLine)^2))    #poverty severity
  F_d[1,7] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^0)))  #poverty incidence
  F_d[1,8] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^1)))  #poverty gap
  F_d[1,9] <- mean(colMeans((exp(y_pred)<PovLine)*(((PovLine-exp(y_pred))/PovLine)^2)))  #poverty severity

  ### Data to return
  return(F_d)
}

NM_Predict <- function(samp, beeta, sig_e, sig_u){

  #Computes input values for NM_Predict
  gamma_nd <- sig_u^2/(sig_e^2+param$nd*sig_u^2)
  alpha    <- sig_e^2
  beta     <- alpha*gamma_nd

  #Computes the predictor for all the small  areas
  y_pred <- t(sapply(1:param$d,FUN=NM_Predict_d,L=param$L,samp=samp,beeta=beeta,sig_e=sig_e,gamma_nd=gamma_nd,alpha=alpha,beta=beta,PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(y_pred)
  names(obj_pred) <- c("I1","G1","S1","I2","G2","S2","I3","G3","S3")
  return(obj_pred)
  
}


### This function predicts the outcome of the non-sampled units

ELL_Population <- function(samp, adjust=TRUE, ELL_method="MOM", resample="TRUE"){

  pop.s  <- as.matrix(samp[samp$Sel%in%1,])
  y_s    <- pop.s[,1]
  x_s    <- pop.s[,3:5]
  lm.s   <- lm(y_s ~ 0 + x_s)
  b_ols  <- as.vector(lm.s$coefficients)
  cov_b  <- vcov(lm.s)
  beeta  <- t(rmvnorm(n=param$L,mean=b_ols,sigma=cov_b))

  if (toupper(resample)=="TRUE") {
    s.ind  <- tapply(1:param$n,INDEX=pop.s[,2],FUN=sample,size=param$nd,replace=TRUE) 
    select <- rep(NA,param$n)
    for (i in 1:param$d){
      select[(1+(i-1)*param$nd):(i*param$nd)] <- as.vector(s.ind[[i]])
    }
    pop.ss <- pop.s[select,]
    y_ss   <- pop.ss[,1]
    x_ss   <- pop.ss[,3:5]
    lm.ss  <- lm(y_ss ~ 0 + x_ss)
  } else if (toupper(resample)=="FALSE") {
      pop.ss <- pop.s
      y_ss   <- y_s
      x_ss   <- x_s
      lm.ss  <- lm.s
    } else {
        print("Error: the parameter resample should be set to TRUE or FALSE") 
        break
    } 

  ## residuals 
  mu_ss      <- (x_ss%*%b_ols)
  res_ols    <- y_ss-mu_ss
  res_ols_d  <- tapply(res_ols, pop.ss[,2],mean)
  res_ols_dj <- res_ols - rep(res_ols_d, each=param$nd)
  
  if (ELL_method=="TRAD") {
    #Area residuals 
    res_d1     <- rep(sample(res_ols_d,size=param$d,replace=TRUE),each=param$Nd)
  
    ## Errors
    res_dj1    <- sample(res_ols_dj,size=param$N,replace=TRUE)
   
    if (adjust == "TRUE"){
      res_d1   <- res_d1 - mean(res_d1)
      res_dj1  <- res_dj1 - rep(tapply(res_dj1, samp[,2],mean),each=param$Nd)
    }
    
    #Computes the predictor for all the small  areas
    y_pred <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d1 + res_dj1  #as.matrix(samp[,3:5])%*%beeta + res_d1 + res_dj1  #traditional method
    
    ### Data to return
    obj_pred        <- data.frame(y_pred)
    names(obj_pred) <- c("pop_ELL_TRAD")
         
    
  } else if (ELL_method=="RES"){
      #Area residuals 
      res_d2     <- rep(res_ols_d,each=param$Nd)
      
      ## Errors
      res_dj2    <- sample(res_ols_dj,size=param$N,replace=TRUE)
         
      if (adjust == "TRUE"){
        res_d2   <- res_d2 - mean(rep(res_ols_d,each=param$Nd)) 
        res_dj2  <- res_dj2 - rep(tapply(res_dj2, samp[,2],mean),each=param$Nd)
      }
      
      #Computes the predictor for all the small  areas
      y_pred <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d2 + res_dj2  #Residual method
      
      ### Data to return
      obj_pred        <- data.frame(y_pred)
      names(obj_pred) <- c("pop_ELL_RES")
      
    } else if (ELL_method=="MOM"){
         y_dev  <- y_ss - rep(tapply(y_ss, pop.ss[,2],mean), each=param$nd)
         x_dev  <- x_ss[,2:3] - cbind(rep(tapply(x_ss[,2],pop.ss[,2],mean),each=param$nd),rep(tapply(x_ss[,3],pop.ss[,2],mean),each=param$nd))
         SSE_1  <- sum(lm(y_dev ~ 0 + x_dev)$residuals^2)
         SSE_2  <- sum(lm.ss$residuals^2)
         den_1  <- nrow(x_ss) - param$d + ncol(x_dev)
         x_mean <- cbind(tapply(x_ss[,1],pop.ss[,2],mean),tapply(x_ss[,2],pop.ss[,2],mean),tapply(x_ss[,3],pop.ss[,2],mean))
         den_2  <- sum(diag(param$nd*(1-param$nd*x_mean%*%solve(t(x_ss)%*%x_ss)%*%t(x_mean))))
       
         sig_em <- max(SSE_1/den_1,1e-9)
         sig_um <- max((SSE_2 - (param$n - 3)*sig_em)/den_2,1e-9)
         gamma_nd <- sig_um^2/(sig_em^2+param$nd*sig_um^2)
         
         u_EB   <- function(d,gamma_nd,res_ols) gamma_nd*sum(res_ols[pop.s[,2]==d])
         
         #Area residuals 
         res_d3     <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$Nd)
         
         res_d3_n   <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$nd)
         res_dj3    <- sample(res_ols-res_d3_n,size=param$N,replace=TRUE)
  
         #Computes the predictor for all the small  areas
         y_pred <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d3 + res_dj3  #MOM method
         
         ### Data to return
	   obj_pred        <- data.frame(y_pred)
         names(obj_pred) <- c("pop_ELL_MOM")
         
      } else {
          print("Error: must indicate the ELL method (TRAD, RES, or MOM)") 
          break
        }

  return(obj_pred)
}


Pov_Predict <- function(y_pred, PovLine){

  ### Complex statistics prediction 
  complex_stats <- function(d, y, PovLine=param$z){
    y_d     <- y[param$areas==d,]
    F_d      <- matrix(data=0, nrow=1, ncol=3)
    F_d[1,1] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^2)))  #poverty severity
    return(F_d)
  }

  F_pred  <- t(sapply(1:param$d,complex_stats, y=y_pred, samp=samp, PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(F_pred)
  names(obj_pred) <- c("I","G","S")
  return(obj_pred)

}

ELL_RES_Predict <- function(samp, adjust=TRUE){
  pop.s  <- as.matrix(samp[samp$Sel%in%1,])
  y_s    <- pop.s[,1]
  x_s    <- pop.s[,3:5]
  lm.s   <- lm(y_s ~ 0 + x_s)
  b_ols  <- as.vector(lm.s$coefficients)
  #cov_b  <- vcov(lm.s)
  #beeta  <- t(rmvnorm(n=param$L,mean=b_ols,sigma=cov_b))

  ## residuals 
  mu_s       <- (x_s%*%b_ols)
  res_ols    <- y_s-mu_s
  res_ols_d  <- tapply(res_ols, pop.s[,2],mean)
  res_ols_dj <- res_ols - rep(res_ols_d, each=param$nd)


  ## Area effects
  res_d      <- rep(res_ols_d,each=param$Nd)

  if (adjust == "TRUE"){
    res_d    <- res_d - mean(rep(res_ols_d,each=param$Nd))        
  }

  ## Errors
  res_dj    <- replicate(param$L,sample(res_ols_dj,size=param$N,replace=TRUE))

  if (adjust == "TRUE"){
    for (l in 1:param$L){
      res_dj[,l]  <- res_dj[,l] - rep(tapply(res_dj[,l], samp[,2],mean),each=param$Nd)
    }
  }

  #Computes the predictor for all the small  areas
  y_pred <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d + res_dj  #Residual method

  ### Complex statistics prediction 
  complex_stats <- function(d, y, samp, PovLine=param$z){
    y_d      <- y[samp[,2]==d,]
    F_d      <- matrix(data=0, nrow=1, ncol=3)
    F_d[1,1] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^2)))  #poverty severity
    return(F_d)
  }

  F_pred  <- t(sapply(1:param$d,complex_stats, y=y_pred, samp=samp, PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(F_pred)
  names(obj_pred) <- c("I","G","S")
  return(obj_pred)

}

ELL_MOM_Predict <- function(samp, adjust=TRUE){
  pop.s  <- as.matrix(samp[samp$Sel%in%1,])
  y_s    <- pop.s[,1]
  x_s    <- pop.s[,3:5]
  lm.s   <- lm(y_s ~ 0 + x_s)
  b_ols  <- as.vector(lm.s$coefficients)
  #cov_b  <- vcov(lm.s)
  #beeta  <- t(rmvnorm(n=param$L,mean=b_ols,sigma=cov_b))

  y_dev  <- y_s - rep(tapply(y_s, pop.s[,2],mean), each=param$nd)
  x_dev  <- x_s[,2:3] - cbind(rep(tapply(x_s[,2],pop.s[,2],mean),each=param$nd),rep(tapply(x_s[,3],pop.s[,2],mean),each=param$nd))
  SSE_1  <- sum(lm(y_dev ~ 0 + x_dev)$residuals^2)
  SSE_2  <- sum(lm.s$residuals^2)
  den_1  <- nrow(x_s) - param$d - ncol(x_dev)
  x_mean <- cbind(tapply(x_s[,1],pop.s[,2],mean),tapply(x_s[,2],pop.s[,2],mean),tapply(x_s[,3],pop.s[,2],mean))
  den_2  <- sum(diag(param$nd*(1-param$nd*x_mean%*%solve(t(x_s)%*%x_s)%*%t(x_mean))))

  sig2_em <- max(SSE_1/den_1,1e-9)
  sig2_um <- max((SSE_2 - (param$n - 3)*sig2_em)/den_2,1e-9)
  gamma_nd <- sig2_um/(sig2_em+param$nd*sig2_um)

  ## residuals 
  mu_s       <- (x_s%*%b_ols)
  res_ols    <- y_s-mu_s
  res_ols_d  <- tapply(res_ols, pop.s[,2],mean)
  res_ols_dj <- res_ols - rep(res_ols_d, each=param$nd)

  u_EB   <- function(d,gamma_nd,res_ols) gamma_nd*sum(res_ols[pop.s[,2]==d])

  ## Area effects
  res_d     <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$Nd)

  if (adjust == "TRUE"){   
    res_d   <- res_d #- mean(res_d)   
  }

  ## Errors
  res_d_n   <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$nd)
  res_dj    <- replicate(param$L,sample(res_ols-res_d_n,size=param$N,replace=TRUE))

  if (adjust == "TRUE"){
    for (l in 1:param$L){
      res_dj[,l]  <- res_dj[,l] #- rep(tapply(res_dj[,l], samp[,2],mean),each=param$Nd)
    }
  }

  #Computes the predictor for all the small  areas
  y_pred <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d + res_dj  #MOM method

  ### Complex statistics prediction 
  complex_stats <- function(d, y, samp, PovLine=param$z){
    y_d     <- y[samp[,2]==d,]
    F_d      <- matrix(data=0, nrow=1, ncol=9)
    F_d[1,1] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_d)<PovLine)*(((PovLine-exp(y_d))/PovLine)^2)))  #poverty severity
    return(F_d)
  }

  F_pred  <- t(sapply(1:param$d,complex_stats, y=y_pred, samp=samp, PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(F_pred)
  names(obj_pred) <- c("I","G","S")
  return(obj_pred)

}


ELL_Predict <- function(samp, adjust=TRUE){
  pop.s  <- as.matrix(samp[samp$Sel%in%1,])
  y_s    <- pop.s[,1]
  x_s    <- pop.s[,3:5]
  lm.s   <- lm(y_s ~ 0 + x_s)
  b_ols  <- as.vector(lm.s$coefficients)
  cov_b  <- vcov(lm.s)
  beeta  <- t(rmvnorm(n=param$L,mean=b_ols,sigma=cov_b))

  y_dev  <- y_s - rep(tapply(y_s, pop.s[,2],mean), each=param$nd)
  x_dev  <- x_s[,2:3] - cbind(rep(tapply(x_s[,2],pop.s[,2],mean),each=param$nd),rep(tapply(x_s[,3],pop.s[,2],mean),each=param$nd))
  SSE_1  <- sum(lm(y_dev ~ 0 + x_dev)$residuals^2)
  SSE_2  <- sum(lm.s$residuals^2)
  den_1  <- nrow(x_s) - param$d - ncol(x_dev)
  x_mean <- cbind(tapply(x_s[,1],pop.s[,2],mean),tapply(x_s[,2],pop.s[,2],mean),tapply(x_s[,3],pop.s[,2],mean))
  den_2  <- sum(diag(param$nd*(1-param$nd*x_mean%*%solve(t(x_s)%*%x_s)%*%t(x_mean))))

  sig2_em <- max(SSE_1/den_1,1e-9)
  sig2_um <- max((SSE_2 - (param$n - 3)*sig2_em)/den_2,1e-9)
  gamma_nd <- sig2_um/(sig2_em+param$nd*sig2_um)

  ## residuals 
  mu_s       <- (x_s%*%b_ols)
  res_ols    <- y_s-mu_s
  res_ols_d  <- tapply(res_ols, pop.s[,2],mean)
  res_ols_dj <- res_ols - rep(res_ols_d, each=param$nd)

  u_EB   <- function(d,gamma_nd,res_ols) gamma_nd*sum(res_ols[pop.s[,2]==d])

  ## Area effects
  res_d1     <- replicate(param$L,rep(sample(res_ols_d,size=param$d,replace=TRUE),each=param$Nd))
  res_d2     <- rep(res_ols_d,each=param$Nd)
  res_d3     <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$Nd)

  if (adjust == "TRUE"){
    res_d1   <- res_d1 - matrix(colMeans(res_d1),param$N,param$L,byrow=TRUE)
    res_d2   <- res_d2 - mean(rep(res_ols_d,each=param$Nd))      
    res_d3   <- res_d3 #- mean(res_d3)   
  }

  ## Errors
  res_dj1    <- replicate(param$L,sample(res_ols_dj,size=param$N,replace=TRUE))
  res_d3_n   <- rep(sapply(1:param$d,u_EB,gamma_nd=gamma_nd,res_ols=res_ols),each=param$nd)
  res_dj3    <- replicate(param$L,sample(res_ols-res_d3_n,size=param$N,replace=TRUE))

  if (adjust == "TRUE"){
    for (l in 1:param$L){
      res_dj1[,l]  <- res_dj1[,l] - rep(tapply(res_dj1[,l], samp[,2],mean),each=param$Nd)
      res_dj3[,l]  <- res_dj3[,l] #- rep(tapply(res_dj3[,l], samp[,2],mean),each=param$Nd)
    }
  }
  res_dj2    <- res_dj1

  #Computes the predictor for all the small  areas
  y_pred1 <- as.matrix(samp[,3:5])%*%beeta            + res_d1 + res_dj1  #traditional method
  y_pred2 <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d2 + res_dj2  #Residual method
  y_pred3 <- as.vector(as.matrix(samp[,3:5])%*%b_ols) + res_d3 + res_dj3  #MOM method

  ### Complex statistics prediction 
  complex_stats <- function(d, y1, y2, y3, samp, PovLine=param$z){
    y_d1     <- y1[samp[,2]==d,]
    y_d2     <- y2[samp[,2]==d,]
    y_d3     <- y3[samp[,2]==d,]
    F_d      <- matrix(data=0, nrow=1, ncol=9)
    F_d[1,1] <- mean(colMeans((exp(y_d1)<PovLine)*(((PovLine-exp(y_d1))/PovLine)^0)))  #poverty incidence
    F_d[1,2] <- mean(colMeans((exp(y_d1)<PovLine)*(((PovLine-exp(y_d1))/PovLine)^1)))  #poverty gap
    F_d[1,3] <- mean(colMeans((exp(y_d1)<PovLine)*(((PovLine-exp(y_d1))/PovLine)^2)))  #poverty severity
    F_d[1,4] <- mean(colMeans((exp(y_d2)<PovLine)*(((PovLine-exp(y_d2))/PovLine)^0)))  #poverty incidence
    F_d[1,5] <- mean(colMeans((exp(y_d2)<PovLine)*(((PovLine-exp(y_d2))/PovLine)^1)))  #poverty gap
    F_d[1,6] <- mean(colMeans((exp(y_d2)<PovLine)*(((PovLine-exp(y_d2))/PovLine)^2)))  #poverty severity
    F_d[1,7] <- mean(colMeans((exp(y_d3)<PovLine)*(((PovLine-exp(y_d3))/PovLine)^0)))  #poverty incidence
    F_d[1,8] <- mean(colMeans((exp(y_d3)<PovLine)*(((PovLine-exp(y_d3))/PovLine)^1)))  #poverty gap
    F_d[1,9] <- mean(colMeans((exp(y_d3)<PovLine)*(((PovLine-exp(y_d3))/PovLine)^2)))  #poverty severity 
    return(F_d)
  }

  F_pred  <- t(sapply(1:param$d,complex_stats, y1=y_pred1, y2=y_pred2, y3=y_pred3, samp=samp, PovLine=param$z))

  ### Data to return
  obj_pred        <- cbind.data.frame(F_pred)
  names(obj_pred) <- c("I1","G1","S1","I2","G2","S2","I3","G3","S3")
  return(obj_pred)

}


