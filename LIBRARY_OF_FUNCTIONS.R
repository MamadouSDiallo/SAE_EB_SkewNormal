############################################################################################################################### 

# THIS PROGRAM CONTAINS FUNCTIONS USED IN THE SKEW-NORMAL SIMULATIONS

# Autheur: Mamadou S. Diallo 
# Date: March 2013

############################################################################################################################### 


### This will load the library where are stored the packages .libPaths('C:/Program
### Files/R/Cont_Library')

### Loading the libraries needed for this program ###

library(Matrix)
library(MASS)
library(mnormt)
library(mvtnorm)
library(numDeriv)
library(sn)


### Remove all objects to start fresh rm(list=ls())

### This function generates a population based on the specs above
Generate_Y <- function(param, adjust = "FALSE") {
    
    ### Model psecifications
    e <- rep(0, param$N)
    
    for (d in 1:param$d) {
        e[(1 + (d - 1) * param$Nd):(d * param$Nd)] <- param$sig_e * (param$delta_e * abs(rnorm(n = param$Nd)) + 
            ((1 - param$delta_e^2)^(1/2)) * rnorm(param$Nd))
    }
    
    ud <- param$sig_u * (param$delta_u * abs(rnorm(n = param$d)) + ((1 - param$delta_u^2)^(1/2)) * 
        rnorm(param$d))
    u <- rep(ud, each = param$Nd)
    
    if (toupper(adjust) == "TRUE") {
        ud <- (ud - param$sig_u * param$delta_u * sqrt(2/pi))/sqrt(1 - 2 * param$delta_u^2/pi)
        u <- rep(ud, each = param$Nd)
        e <- (e - param$sig_e * param$delta_e * sqrt(2/pi))/sqrt(1 - 2 * param$delta_e^2/pi)
        # u <- u - mean(u) e <- e - mean(e)
    }
    
    ### Outcome variable Y
    Y <- param$X %*% param$beeta + u + e
    
    ### Data to return
    obj <- cbind.data.frame(Y, param$areas, param$X)
    names(obj) = c("Y", "Area", "X0", "X1", "X2")
    return(list(dat = obj, u = ud, e = e))
}


### This function selects the sample in each areas independently
samp_Y <- function(Y, param, sys = FALSE, proportional = FALSE) {
    # set.seed(param$d*param$Nd*param$nd*i)
    samp <- matrix(0, nrow = param$n, ncol = 1)
    sel <- rep(0, param$N)
    weight <- rep(0, param$N)
    if (sys) {
        for (d in 1:param$d) {
            range <- param$Nd/param$nd
            starting <- runif(1, min = 1 + (d - 1) * param$Nd, max = (d - 1) * param$Nd + range)
            samp[(1 + (d - 1) * param$nd):(d * param$nd), 1] <- round(seq(from = starting, to = d * 
                param$Nd, by = range), 0)
        }
        Y_sort <- Y[order(Y$Area, Y$Y), ]
        sel[samp[1:param$n, 1]] <- 1
        weight[samp[1:param$n, 1]] <- rep(param$N/param$n, param$n)
        obj <- cbind.data.frame(Y_sort, sel, weight)
    } else {
        for (d in 1:param$d) {
            if (proportional) {
                prob_select = Y$Y[(1 + (d - 1) * param$Nd):(d * param$Nd)]
            } else {
                prob_select = rep(1, param$Nd)
            }
            samp[(1 + (d - 1) * param$nd):(d * param$nd), 1] <- sample((1:param$N)[Y[, 2] == d], 
                param$nd, prob = prob_select)
        }
        sel[samp[1:param$n, 1]] <- 1
        weight[samp[1:param$n, 1]] <- rep(param$N/param$n, param$n)
        obj <- cbind.data.frame(Y, sel, weight)
    }
    names(obj) = c("Y", "Area", "X0", "X1", "X2", "Sel", "Weight")
    return(obj)
}

### This function shows the empirical selection bias
selection_bias <- function(param, adjust = "TRUE", sim) {
    bias <- matrix(NA, sim, 3)
    for (i in 1:sim) {
        y <- Generate_Y(param, adjust = adjust)
        s1 <- samp_Y(y, param, sys = FALSE, proportional = FALSE)
        s2 <- samp_Y(y, param, sys = TRUE, proportional = FALSE)
        s3 <- samp_Y(y, param, sys = FALSE, proportional = TRUE)
        bias[i, 1] <- 100 * (mean(s1$Y[s1$Sel == 1]) - mean(y$Y))/mean(y$Y)
        bias[i, 2] <- 100 * (mean(s2$Y[s2$Sel == 1]) - mean(y$Y))/mean(y$Y)
        bias[i, 3] <- 100 * (mean(s3$Y[s3$Sel == 1]) - mean(y$Y))/mean(y$Y)
    }
    obj <- data.frame(bias)
    names(obj) <- c("SRS", "SYSTEMATIC", "PROPORTIONAL")
    
    return(obj)
}

# colMeans(selection_bias(param,adjust='TRUE', sim=100))

### This function returns the parameters of the distribution of Yd
joint_dist <- function(mu_u = 0, sig_u, lambda_u, mu_e = 0, sig_e, lambda_e) {
    mu_u <- mu_u
    mu_e <- mu_e
    mu_yd <- (mu_u + mu_e) * rep(1, param$nd)
    C_yd <- 1/((pnorm(0, 0, 1 + lambda_e^2)^(param$nd)) * pnorm(0, 0, 1 + lambda_u^2))
    Sig_u <- (sig_u^2) * rep(1, param$nd) %*% t(rep(1, param$nd))
    Sig_e <- (sig_e^2) * diag(param$nd)
    Sig_yd <- Sig_e + Sig_u
    Sig_yd_inv <- (1/sig_e^2) * (diag(param$nd) - (sig_u^2/(sig_e^2 + param$nd * sig_u^2)) * rep(1, 
        param$nd) %*% t(rep(1, param$nd)))
    D_u <- (lambda_u/(param$nd * sig_u)) * t(rep(1, param$nd))
    D_e <- (lambda_e/sig_e) * diag(param$nd)
    D_yd <- t(cbind(Sig_e %*% t(D_e), Sig_u %*% t(D_u))) %*% Sig_yd_inv
    nu_yd <- rep(0, param$nd + 1)
    Gamma_yd <- as.matrix(diag(param$nd + 1) + bdiag(D_e, D_u) %*% bdiag(Sig_e, Sig_u) %*% t(bdiag(D_e, 
        D_u)) - rbind(D_e %*% Sig_e, D_u %*% Sig_u) %*% Sig_yd_inv %*% t(rbind(D_e %*% Sig_e, D_u %*% 
        Sig_u)))
    return(list(C_Yd = C_yd, mu_Yd = mu_yd, Sigma_Yd = Sig_yd, Sigma_Yd_inv = Sig_yd_inv, D_Yd = D_yd, 
        nu_Yd = nu_yd, Gamma_Yd = Gamma_yd))
}

### This function returns the parameters of the distribution of Ydr given Yds
cond_dist <- function(mu_u = 0, sig_u, lambda_u, mu_e = 0, sig_e, lambda_e, XdrBeta, RESds) {
    mu_u <- mu_u
    mu_e <- mu_e
    Ndr <- param$Nd - param$nd
    Nds <- param$nd
    Nd <- Ndr + Nds
    mu_ydrs <- XdrBeta + (mu_u + mu_e) * rep(1, Ndr) + (sig_u^2/(sig_e^2 + Nds * sig_u^2)) * (rep(1, 
        Nds) %*% RESds)[1] * rep(1, Ndr)
    Sig_ydrs <- (sig_e^2) * (diag(Ndr) + (sig_u^2/(sig_e^2 + Nds * sig_u^2)) * rep(1, Ndr) %*% t(rep(1, 
        Ndr)))
    D_u <- matrix((lambda_u/sig_u) * (sig_u^2/(sig_e^2 + Nd * sig_u^2)), 1, Ndr)
    D_e <- (lambda_e/sig_e) * (rbind(diag(Ndr), matrix(0, Nds, Ndr)) - (sig_u^2/(sig_e^2 + Nd * sig_u^2)) * 
        rep(1, Nd) %*% t(rep(1, Ndr)))
    D_ydrs <- as.matrix(rbind(D_e, D_u))
    nu_ydrs <- -as.matrix(rbind(matrix(0, Ndr, Nds), (lambda_e/sig_e) * (diag(Nds) - (sig_u^2/(sig_e^2 + 
        Nds * sig_u^2)) * rep(1, Nds) %*% t(rep(1, Nds))), matrix((lambda_u/sig_u) * (sig_u^2/(sig_e^2 + 
        Nds * sig_u^2)), 1, Nds)) %*% RESds)
    Gamma_ydrs <- joint_dist(mu_u = mu_u, sig_u = sig_u, lambda_u = lambda_u, mu_e = mu_e, sig_e = sig_e, 
        lambda_e = lambda_e)$Gamma_Yd
    
    return(list(mu_Ydrs = mu_ydrs, Sigma_Ydrs = Sig_ydrs, D_Ydrs = D_ydrs, nu_Ydrs = nu_ydrs, Gamma_Ydrs = Gamma_ydrs))
}


### This function computes an approximation of Phi for small are d Reference: volume1: Models and
### Applications 2nd ed. author: Samuel Kotz, N. Balakrishnan and N. L. Jonhson

# approx.Phi.d <- function(n.int=1000,ll=-10,uu=10,Joint.Dist,d,beeta){ I <- diag(param$nd) J <-
# matrix(1,param$nd,param$nd) Area <- samp$Area[samp$Sel==1] y <- samp$Y[samp$Sel==1] x <-
# cbind(samp$X0,samp$X1,samp$X2)[samp$Sel==1,] mu <- x%*%beeta+Joint.Dist$mu_Yd R <-
# diag(1/sqrt(diag(Joint.Dist$Gamma_Yd)))%*%Joint.Dist$Gamma_Yd%*%diag(1/sqrt(diag(Joint.Dist$Gamma_Yd)))
# cj <- as.vector(c(rep(sqrt(R[1,param$nd]),param$nd),R[1,param$nd+1]/sqrt(R[1,param$nd]))) zj <-
# as.vector((1/sqrt(diag(Joint.Dist$Gamma_Yd)))*(Joint.Dist$D_Yd%*%((y-mu)[Area==d]))-Joint.Dist$nu_Yd)
# u0 <- seq(ll,uu*((n.int-1)/n.int),length.out=n.int+1)+(uu-ll)/(2*n.int) dim.u0 <- length(u0)
# phi.hj <- function(u0,zj,cj,k) prod(pnorm((zj-cj*u0[k])/sqrt(1-cj^2))) phi <-
# sum(dnorm(u0)*sapply(1:dim.u0,FUN=phi.hj,u0=u0,zj=zj,cj=cj))*((uu-ll)/(n.int+1)) return(phi) }

# approx.Phi.d <- function(u0,Joint.Dist,d,beeta){ I <- diag(param$nd) J <-
# matrix(1,param$nd,param$nd) Area <- samp$Area[samp$Sel==1] y <- samp$Y[samp$Sel==1] x <-
# cbind(samp$X0,samp$X1,samp$X2)[samp$Sel==1,] mu <- x%*%beeta+mu_Yd R <-
# diag(1/sqrt(diag(Joint.Dist$Gamma_Yd)))%*%Joint.Dist$Gamma_Yd%*%diag(1/sqrt(diag(Joint.Dist$Gamma_Yd)))
# cj <- as.vector(c(rep(sqrt(R[1,param$nd]),param$nd),R[1,param$nd+1]/sqrt(R[1,param$nd]))) zj <-
# as.vector((1/sqrt(diag(Joint.Dist$Gamma_Yd)))*(Joint.Dist$D_Yd%*%((y-mu)[Area==d]))-Joint.Dist$nu_Yd)
# phi.j <- function(u0,zj,cj) dnorm(u0)*prod(pnorm((zj-cj*u0)/sqrt(1-cj^2))) return(mapply(phi.j,
# u0, MoreArgs=list(zj=zj,cj=cj), SIMPLIFY = TRUE)) }


### This function computes the log-likelihoog for the LMM
log.like <- function(par, adjust = "FALSE", model = "ue") {
    
    if (par[4] <= 0 | par[5] <= 0) {
        loglik <- -1e+06
    } else {
        beeta <- par[1:3]
        sig2_u <- par[4]
        sig2_e <- par[5]
        ## given the specified model, lambda_u or lambda_e can be zero
        if (model == "u") {
            lambda_u <- par[6]
            lambda_e <- 0
        } else if (model == "e") {
            lambda_u <- 0
            lambda_e <- par[6]
        } else if (model == "ue") {
            lambda_u <- par[6]
            lambda_e <- par[7]
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
        delta_u <- lambda_u/sqrt(1 + lambda_u^2)
        delta_e <- lambda_e/sqrt(1 + lambda_e^2)
        Area <- samp$Area[samp$Sel == 1]
        y <- samp$Y[samp$Sel == 1]
        x <- cbind(samp$X0, samp$X1, samp$X2)[samp$Sel == 1, ]
        if (adjust == "TRUE") {
            mu_u <- -delta_u * sqrt(sig2_u * 2/pi)
            mu_e <- -delta_e * sqrt(sig2_e * 2/pi)
        } else {
            mu_u <- mu_e <- 0
        }
        joint <- joint_dist(mu_u = mu_u, mu_e = mu_e, sig_u = sqrt(sig2_u), lambda_u = lambda_u, 
            sig_e = sqrt(sig2_e), lambda_e = lambda_e)
        mu <- x %*% beeta + rep(joint$mu_Yd, param$d)
        R <- diag(1/sqrt(diag(joint$Gamma_Yd))) %*% joint$Gamma_Yd %*% diag(1/sqrt(diag(joint$Gamma_Yd)))
        if (model == "u") {
            cj <- rep(0, param$nd + 1)
        } else if (model == "e") {
            cj <- as.vector(c(rep(sqrt(R[1, param$nd]), param$nd), 0))
        } else if (model == "ue") {
            cj <- as.vector(c(rep(sqrt(R[1, param$nd]), param$nd), R[1, param$nd + 1]/sqrt(R[1, param$nd])))
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
        
        
        ### The function Phi.d computes an approximation of Phi for small are d Reference: volume1: Models
        ### and Applications 2nd ed. author: Samuel Kotz, N. Balakrishnan and N. L. Jonhson
        Phi.d <- function(u0, zj, cj) {
            phi.j <- function(u0, zj, cj) dnorm(u0) * prod(pnorm((zj - cj * u0)/sqrt(1 - cj^2)))
            return(mapply(phi.j, u0, MoreArgs = list(zj = zj, cj = cj), SIMPLIFY = TRUE))
        }
        
        ## This functon computes the multivariate probability Phi
        log.like.d <- function(d, joint, lambda_u, lambda_e, y, mu, Area, cj) {
            # limit and nd.points were chosen to give an excellent approximation more than 20 digits after
            # the dot
            limit <- 25
            nb.points <- 1000
            int.length <- 2 * limit/nb.points
            start.int <- -limit + int.length * 0.5
            end.int <- limit - int.length * 0.5
            univers <- seq(from = start.int, to = end.int, by = int.length)
            
            zj <- as.vector((1/sqrt(diag(joint$Gamma_Yd))) * (joint$D_Yd %*% ((y - mu)[Area == d])) - 
                joint$nu_Yd)
            logPhi <- loglik <- numeric(1)
            
            if (lambda_e == 0 & lambda_u == 0) {
                logPhi <- -(param$nd + 1) * log(2)
                # } else if (lambda_e == 0){ logPhi <-
                # -param$nd*log(2)+pnorm(as.numeric(joint$D_Yd%*%((y-mu)[Area==d]))[param$nd+1],
                # sd=sqrt(joint$Gamma_Yd[param$nd+1,param$nd+1]),log=TRUE)
            } else {
                logPhi <- log(int.length * sum(Phi.d(u0 = univers, zj = zj, cj = cj)))
            }
            loglik <- -0.5 * (determinant(joint$Sigma_Yd)$modulus[1] + t((y - mu)[Area == d]) %*% 
                solve(joint$Sigma_Yd, (y - mu)[Area == d])) + logPhi
            return(loglik)
        }
        ## Sum the area's log likelihood to get the overall likelihood
        loglik <- sum(sapply(1:param$d, FUN = log.like.d, joint = joint, lambda_u = lambda_u, lambda_e = lambda_e, 
            y = y, mu = mu, Area = Area, cj = cj))
    }
    return(-loglik[1])
}


### This function computes components of the derivarive of Phi
dev.Phi.d <- function(Joint.Dist, d, beeta, sigma2_u, sigma2_e, lambda_u, lambda_e, delta_u, delta_e, 
    Area, y, x, mu, R, cj, adjust, model) {
    
    limit <- 10
    nb.points <- 100  # this is a good enough approximation for the derivatives
    int.length <- 2 * limit/nb.points
    start.int <- -limit + int.length * 0.5
    end.int <- limit - int.length * 0.5
    univers <- seq(from = start.int, to = end.int, by = int.length)
    
    zj <- as.vector(diag(1/sqrt(diag(Joint.Dist$Gamma_Yd))) %*% Joint.Dist$D_Yd %*% (y - mu)[Area == 
        d] - Joint.Dist$nu_Yd)
    mu.beta <- x
    
    if (adjust == "TRUE") {
        mu.sigma2.u <- -rep(delta_u * sqrt(1/(2 * pi * sigma2_u)), param$nd)
        mu.sigma2.e <- -rep(delta_e * sqrt(1/(2 * pi * sigma2_e)), param$nd)
        if (model == "u") {
            mu.lambda.u <- -rep(((1 + lambda_u^2)^(-3/2)) * sqrt(2 * sigma2_u/pi), param$nd)
        } else if (model == "e") {
            mu.lambda.e <- -rep(((1 + lambda_e^2)^(-3/2)) * sqrt(2 * sigma2_e/pi), param$nd)
        } else if (model == "ue") {
            mu.lambda.u <- -rep(((1 + lambda_u^2)^(-3/2)) * sqrt(2 * sigma2_u/pi), param$nd)
            mu.lambda.e <- -rep(((1 + lambda_e^2)^(-3/2)) * sqrt(2 * sigma2_e/pi), param$nd)
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
    } else {
        mu.sigma2.u <- rep(0, param$nd)
        mu.sigma2.e <- rep(0, param$nd)
        if (model == "u") {
            mu.lambda.u <- rep(0, param$nd)
        } else if (model == "e") {
            mu.lambda.e <- rep(0, param$nd)
        } else if (model == "ue") {
            mu.lambda.u <- rep(0, param$nd)
            mu.lambda.e <- rep(0, param$nd)
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
    }
    sigma.u.e <- sigma2_e + param$nd * sigma2_u
    sigma.delta.1 <- sigma2_e + param$nd * sigma2_u + (lambda_e^2) * sigma2_u
    sigma.delta.2 <- sigma2_e + param$nd * sigma2_u + (lambda_u^2) * sigma2_e
    c.j.beta <- matrix(0, param$nd + 1, length(beeta))
    c.j.sigma2.u <- c(rep(0.5 * (lambda_e/sqrt(sigma2_u)) * sigma2_e * (sigma.delta.1^(-3/2)), param$nd), 
        0.5 * (lambda_u/sqrt(sigma2_e)) * param$nd * sigma2_e * (sigma.delta.2^(-3/2)))
    c.j.sigma2.e <- c(rep(-0.5 * (lambda_e * sqrt(sigma2_u)) * (sigma.delta.1^(-3/2)), param$nd), 
        -0.5 * (lambda_u/sqrt(sigma2_e)) * param$nd * sigma2_u * (sigma.delta.2^(-3/2)))
    gamma.j <- c(rep(sqrt((sigma.u.e)/sigma.delta.1), param$nd), sqrt((sigma.u.e)/sigma.delta.2))
    gamma.j.beta <- rep(0, param$nd + 1)
    gamma.j.sigma2.u <- c(rep(-0.5 * (1/sqrt(sigma.u.e)) * ((lambda_e^2) * sigma2_e)/(sigma.delta.1^(3/2)), 
        param$nd), 0.5 * (1/sqrt(sigma.u.e)) * (param$nd * (lambda_u^2) * sigma2_e)/(sigma.delta.2^(3/2)))
    gamma.j.sigma2.e <- c(rep(0.5 * (1/sqrt(sigma.u.e)) * ((lambda_e^2) * sigma2_u)/(sigma.delta.1^(3/2)), 
        param$nd), -0.5 * (1/sqrt(sigma.u.e)) * (param$nd * (lambda_u^2) * sigma2_u)/(sigma.delta.2^(3/2)))
    D.j.beta <- matrix(0, nrow = param$nd + 1, ncol = param$nd)
    D.j.sigma2.u <- rbind(-lambda_e * sqrt(sigma2_e) * ((sigma.u.e)^(-2)) * rep(1, param$nd) %*% 
        t(rep(1, param$nd)), 0.5 * (lambda_u/sqrt(sigma2_u)) * ((sigma2_e - param$nd * sigma2_u) * 
        ((sigma.u.e)^(-2))) * t(rep(1, param$nd)))
    D.j.sigma2.e <- rbind(-0.5 * (lambda_e/(sigma2_e^(3/2))) * (diag(param$nd) - (sigma2_u * (3 * 
        sigma2_e + param$nd * sigma2_u)/(sigma.u.e^2)) * rep(1, param$nd) %*% t(rep(1, param$nd))), 
        -(lambda_u * sqrt(sigma2_u)/(sigma.u.e^2)) * t(rep(1, param$nd)))
    z.j.beta <- gamma.j * (Joint.Dist$D_Yd %*% (-mu.beta[Area == d, ]))
    z.j.sigma2.u <- as.vector(gamma.j.sigma2.u * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + gamma.j * 
        (D.j.sigma2.u %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.sigma2.u)))
    z.j.sigma2.e <- as.vector(gamma.j.sigma2.e * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + gamma.j * 
        (D.j.sigma2.e %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.sigma2.e)))
    
    if (model == "u") {
        c.j.lambda.u <- c(rep(0, param$nd), -sqrt(sigma2_e) * sigma.u.e * (sigma.delta.2^(-3/2)))
        gamma.j.lambda.u <- c(rep(0, param$nd), -lambda_u * sigma2_e * sqrt(sigma.u.e) * (sigma.delta.2^(-3/2)))
        D.j.lambda.u <- rbind(0 * diag(param$nd), (sqrt(sigma2_u)/(sigma.u.e)) * t(rep(1, param$nd)))
        z.j.lambda.u <- as.vector(gamma.j.lambda.u * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + 
            gamma.j * (D.j.lambda.u %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.lambda.u)))
    } else if (model == "e") {
        c.j.lambda.e <- c(rep(sqrt(sigma2_u) * sigma.u.e * (sigma.delta.1^(-3/2)), param$nd), 0)
        gamma.j.lambda.e <- c(rep(-lambda_e * sigma2_u * sqrt(sigma.u.e) * (sigma.delta.1^(-3/2)), 
            param$nd), 0)
        D.j.lambda.e <- rbind((1/sqrt(sigma2_e)) * (diag(param$nd) - (sigma2_u/(sigma.u.e)) * rep(1, 
            param$nd) %*% t(rep(1, param$nd))), 0 * t(rep(1, param$nd)))
        z.j.lambda.e <- as.vector(gamma.j.lambda.e * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + 
            gamma.j * (D.j.lambda.e %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.lambda.e)))
    } else if (model == "ue") {
        c.j.lambda.u <- c(rep(0, param$nd), -sqrt(sigma2_e) * sigma.u.e * (sigma.delta.2^(-3/2)))
        c.j.lambda.e <- c(rep(sqrt(sigma2_u) * sigma.u.e * (sigma.delta.1^(-3/2)), param$nd), 0)
        gamma.j.lambda.u <- c(rep(0, param$nd), -lambda_u * sigma2_e * sqrt(sigma.u.e) * (sigma.delta.2^(-3/2)))
        gamma.j.lambda.e <- c(rep(-lambda_e * sigma2_u * sqrt(sigma.u.e) * (sigma.delta.1^(-3/2)), 
            param$nd), 0)
        D.j.lambda.u <- rbind(0 * diag(param$nd), (sqrt(sigma2_u)/(sigma.u.e)) * t(rep(1, param$nd)))
        D.j.lambda.e <- rbind((1/sqrt(sigma2_e)) * (diag(param$nd) - (sigma2_u/(sigma.u.e)) * rep(1, 
            param$nd) %*% t(rep(1, param$nd))), 0 * t(rep(1, param$nd)))
        z.j.lambda.u <- as.vector(gamma.j.lambda.u * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + 
            gamma.j * (D.j.lambda.u %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.lambda.u)))
        z.j.lambda.e <- as.vector(gamma.j.lambda.e * (Joint.Dist$D_Yd %*% ((y - mu)[Area == d])) + 
            gamma.j * (D.j.lambda.e %*% ((y - mu)[Area == d])) + gamma.j * (Joint.Dist$D_Yd %*% (-mu.lambda.e)))
    } else {
        print("Please specify a valid model (u, e, or ue)")
        break
    }
    
    ### The function Phi.d computes an approximation of Phi for small are d Reference: volume1: Models
    ### and Applications 2nd ed. author: Samuel Kotz, N. Balakrishnan and N. L. Jonhson
    Phi.d <- function(u0, zj, cj) {
        phi.j <- function(u0, zj, cj) dnorm(u0) * prod(pnorm((zj - cj * u0)/sqrt(1 - cj^2)))
        return(mapply(phi.j, u0, MoreArgs = list(zj = zj, cj = cj), SIMPLIFY = TRUE))
    }
    
    Phi <- int.length * sum(Phi.d(u0 = univers, zj = zj, cj = cj))
    
    c.j.coef <- 1/sqrt(1 - cj^2)
    
    ### This is the derivative of Phi.d
    Phi.dd <- function(u0, coef, zj, cj, cj.dev, zj.dev) {
        phi.jj <- function(u0, coef, zj, cj, zj.dev, cj.dev) {
            hj <- (zj - cj * u0)/sqrt(1 - cj^2)
            Phi <- pnorm(hj)
            small <- Phi < 9.98012604599318e-322
            Phi[small] <- 9.98012604599318e-322
            result <- dnorm(u0) * sum(coef * ((coef^2) * (zj * cj - u0) * cj.dev + zj.dev) * dnorm(hj) * 
                as.vector(prod(Phi, na.rm = TRUE)/Phi))
            return(result)
        }
        return(mapply(phi.jj, u0, MoreArgs = list(coef = coef, zj = zj, cj = cj, zj.dev = zj.dev, 
            cj.dev = cj.dev), SIMPLIFY = TRUE))
    }
    
    dev.coef1.beta <- as.vector(t(x[Area == d, ]) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - mu)[Area == 
        d]))
    
    if (adjust == "TRUE") {
        dev.coef1.sigma2.u <- -2 * as.vector(t(mu.sigma2.u) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
            mu)[Area == d]))
        dev.coef1.sigma2.e <- -2 * as.vector(t(mu.sigma2.e) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
            mu)[Area == d]))
        if (model == "u") {
            dev.coef1.lambda.u <- -2 * as.vector(t(mu.lambda.u) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
                mu)[Area == d]))
        } else if (model == "e") {
            dev.coef1.lambda.e <- -2 * as.vector(t(mu.lambda.e) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
                mu)[Area == d]))
        } else if (model == "ue") {
            dev.coef1.lambda.e <- -2 * as.vector(t(mu.lambda.e) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
                mu)[Area == d]))
            dev.coef1.lambda.u <- -2 * as.vector(t(mu.lambda.u) %*% Joint.Dist$Sigma_Yd_inv %*% ((y - 
                mu)[Area == d]))
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
    } else {
        dev.coef1.sigma2.u <- 0
        dev.coef1.sigma2.e <- 0
        if (model == "u") {
            dev.coef1.lambda.u <- 0
        } else if (model == "e") {
            dev.coef1.lambda.e <- 0
        } else if (model == "ue") {
            dev.coef1.lambda.e <- 0
            dev.coef1.lambda.u <- 0
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
    }
    
    dev.coef2.u <- as.double(-0.5 * (sum(diag(Joint.Dist$Sigma_Yd_inv %*% (rep(1, param$nd) %*% t(rep(1, 
        param$nd))))) - t((y - mu)[Area == d]) %*% Joint.Dist$Sigma_Yd_inv %*% (rep(1, param$nd) %*% 
        t(rep(1, param$nd))) %*% Joint.Dist$Sigma_Yd_inv %*% as.vector((y - mu)[Area == d])))
    dev.coef2.e <- as.double(-0.5 * (sum(diag(Joint.Dist$Sigma_Yd_inv)) - t((y - mu)[Area == d]) %*% 
        Joint.Dist$Sigma_Yd_inv %*% Joint.Dist$Sigma_Yd_inv %*% ((y - mu)[Area == d])))
    
    phi.j.beta0 <- dev.coef1.beta[1] + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, zj = zj, 
        cj = cj, cj.dev = rep(0, param$nd + 1), zj.dev = as.vector(z.j.beta[, 1])))/Phi
    phi.j.beta1 <- dev.coef1.beta[2] + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, zj = zj, 
        cj = cj, cj.dev = rep(0, param$nd + 1), zj.dev = as.vector(z.j.beta[, 2])))/Phi
    phi.j.beta2 <- dev.coef1.beta[3] + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, zj = zj, 
        cj = cj, cj.dev = rep(0, param$nd + 1), zj.dev = as.vector(z.j.beta[, 3])))/Phi
    
    phi.j.sigma2.u <- -0.5 * dev.coef1.sigma2.u + dev.coef2.u + int.length * sum(Phi.dd(u0 = univers, 
        coef = c.j.coef, zj = zj, cj = cj, cj.dev = c.j.sigma2.u, zj.dev = z.j.sigma2.u))/Phi
    phi.j.sigma2.e <- -0.5 * dev.coef1.sigma2.e + dev.coef2.e + int.length * sum(Phi.dd(u0 = univers, 
        coef = c.j.coef, zj = zj, cj = cj, cj.dev = c.j.sigma2.e, zj.dev = z.j.sigma2.e))/Phi
    
    if (model == "u") {
        phi.j.lambda.u <- -0.5 * dev.coef1.lambda.u + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, 
            zj = zj, cj = cj, cj.dev = c.j.lambda.u, zj.dev = z.j.lambda.u))/Phi
        gradient <- c(phi.j.beta0, phi.j.beta1, phi.j.beta2, phi.j.sigma2.u, phi.j.sigma2.e, phi.j.lambda.u)
    } else if (model == "e") {
        phi.j.lambda.e <- -0.5 * dev.coef1.lambda.e + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, 
            zj = zj, cj = cj, cj.dev = c.j.lambda.e, zj.dev = z.j.lambda.e))/Phi
        gradient <- c(phi.j.beta0, phi.j.beta1, phi.j.beta2, phi.j.sigma2.u, phi.j.sigma2.e, phi.j.lambda.e)
    } else if (model == "ue") {
        phi.j.lambda.u <- -0.5 * dev.coef1.lambda.u + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, 
            zj = zj, cj = cj, cj.dev = c.j.lambda.u, zj.dev = z.j.lambda.u))/Phi
        phi.j.lambda.e <- -0.5 * dev.coef1.lambda.e + int.length * sum(Phi.dd(u0 = univers, coef = c.j.coef, 
            zj = zj, cj = cj, cj.dev = c.j.lambda.e, zj.dev = z.j.lambda.e))/Phi
        gradient <- c(phi.j.beta0, phi.j.beta1, phi.j.beta2, phi.j.sigma2.u, phi.j.sigma2.e, phi.j.lambda.u, 
            phi.j.lambda.e)
    } else {
        print("Please specify a valid model (u, e, or ue)")
        break
    }
    
    return(gradient)
}


### This function provides the analytical gradient
grad.log.like <- function(par, adjust = "FALSE", model = "ue") {
    if (par[4] <= 0 | par[5] <= 0) {
        dev <- rep(-1e+06, length(par))
    } else {
        beeta <- par[1:3]
        sig2_u <- par[4]
        sig2_e <- par[5]
        ## given the specified model, lambda_u or lambda_e can be zero
        if (model == "u") {
            lambda_u <- par[6]
            lambda_e <- 0
        } else if (model == "e") {
            lambda_u <- 0
            lambda_e <- par[6]
        } else if (model == "ue") {
            lambda_u <- par[6]
            lambda_e <- par[7]
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
        delta_u <- lambda_u/sqrt(1 + lambda_u^2)
        delta_e <- lambda_e/sqrt(1 + lambda_e^2)
        Area <- samp$Area[samp$Sel == 1]
        y <- samp$Y[samp$Sel == 1]
        x <- cbind(samp$X0, samp$X1, samp$X2)[samp$Sel == 1, ]
        
        if (adjust == "TRUE") {
            mu_u <- -delta_u * sqrt(sig2_u * 2/pi)
            mu_e <- -delta_e * sqrt(sig2_e * 2/pi)
        } else {
            mu_u <- mu_e <- 0
        }
        
        joint <- joint_dist(mu_u = mu_u, mu_e = mu_e, sig_u = sqrt(sig2_u), lambda_u = lambda_u, 
            sig_e = sqrt(sig2_e), lambda_e = lambda_e)
        
        mu <- x %*% beeta + joint$mu_Yd
        dev <- rep(0, length(par))
        R <- diag(1/sqrt(diag(joint$Gamma_Yd))) %*% joint$Gamma_Yd %*% diag(1/sqrt(diag(joint$Gamma_Yd)))
        if (model == "u") {
            cj <- rep(0, param$nd + 1)
        } else if (model == "e") {
            cj <- as.vector(c(rep(sqrt(R[1, param$nd]), param$nd), 0))
        } else if (model == "ue") {
            cj <- as.vector(c(rep(sqrt(R[1, param$nd]), param$nd), R[1, param$nd + 1]/sqrt(R[1, param$nd])))
        } else {
            print("Please specify a valid model (u, e, or ue)")
            break
        }
        
        dev <- rowSums(sapply(1:param$d, FUN = dev.Phi.d, Joint.Dist = joint, beeta = beeta, sigma2_u = sig2_u, 
            sigma2_e = sig2_e, lambda_u = lambda_u, lambda_e = lambda_e, delta_u = delta_u, delta_e = delta_e, 
            Area = Area, y = y, x = x, mu = mu, R = R, cj = cj, adjust = adjust, model = model), 
            na.rm = TRUE)
    }
    
    return(-dev)
}


### This function provides a numerical proximation of the analytical gradient
approx.grad.log.like <- function(par, method = "simple", adjust = "TRUE", model = "ue") {
    
    return(grad(func = log.like, x = par, method = method, adjust = adjust, model = model))
    
}






