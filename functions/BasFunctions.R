# Date: 2021/03/08
####Semiparametric efficient estimation of non-ignorable missing data
#########################################################################
##
## R code to implement the Semi-parametric efficient estimator of missing not at random with shadow variables 
## Goal: Study the performance of the proposed methods on estimation of the outcome mean
##
## Methods: Inverse Probability weighting;
##          Doubly robust estimators [2016, Miao et al.];
##          Semi-parametric efficient estimator
##
## Notation:
##  Y = outcome with missing not at random
##  X = fully observed covariates
##    Xs = covariates used in propensity score
##    Xr = covariates used in outcome model for Y
##    xz = covariates used in outcome model for Z
##  Z = shadow variables: independent with missingness indicator given Y and X
##  R = Missingness indicator

library(numDeriv)

#########################################################################
### Function: Generalized Method Moments Estimation
GMM <- function(f, para, data){
    GMM.colmean <- apply(f(para=para, data=data), 2, mean)
    GMM.Func <- sum(GMM.colmean^2)  
    
    return(GMM.Func)
}

### Inverse probability weight
INVweight <- function(Xs, Y, R, gamma, alpha){
    invw <- numeric(length(Y))
    alpha0 <- cbind(Y,Xs) %*% c(gamma, alpha)
    invw[R==1] <- 1/plogis(alpha0[R==1])
    return(invw)
} 

#####################################################################
### Estimating equations for inverse probability weighting
# Input: para = (psi, gamma, alpha)
# Output: A n*(ds+2) matrix, each column stands for one estimating equation
IPWmrf <-  function(para, data){
    Y <- data$Y
    R <- data$R
    Z <- data$Z
    Xs <- data$Xs
    psi <- para[1]
    gamma <- para[2]
    alpha <- para[-c(1,2)] 
    
    w <- INVweight(Xs,Y,R,gamma,alpha)
    g1 <- (w * R - 1) * cbind(Xs, Z)  # Choose user-specified function h(X,Z) exponential
    g2 <- rep(-psi, length(Y))
    g2[R==1] <- -psi + w[R==1] * Y[R==1] 
    g <- cbind(g1,g2)
    
    return(g)
}
#####################################################################
### Estimating equations for regression based estimation
# Input: para = (psi, gamma, beta1, beta2, sigmay2, sigmaz2)
# Output: A n*(dr+dz+5) matrix, each column stands for one estimating equation
REGmrf <-  function(para, data){
    Y <- data$Y
    R <- data$R
    Y[R==0] <- 0
    Z <- data$Z 
    Xr <- data$Xr
    Xz <- data$Xz
    dr <- dim(Xr)[2]
    dz <- dim(Xz)[2]
    
    psi <- para[1]
    gamma <- para[2]
    beta1 <- para[3:(dr+3)]
    beta2 <- para[(dr+4):(dr+dz+3)]
    beta12 <- beta1[1]
    sigmay2 <- para[dr+dz+4]
    sigmaz2 <- para[dr+dz+5]
    
    Sy <- cbind(as.numeric(Y - cbind(Z,Xr) %*% beta1) * cbind(Z,Xr),  
                (Y - cbind(Z,Xr) %*% beta1)^2 - sigmay2) / mean(R)
    Sz <- cbind(as.numeric(Z - Xz %*% beta2)*Xz, 
                (Z - Xz %*% beta2)^2 - sigmaz2) / mean(R)
    Sy[R==0,] <- 0
    Sz[R==0,] <- 0
    
    beta0 <- Xz %*% beta2 - gamma * beta12 * sigmaz2
    g1 <- (1-R) * (Z - beta0)
    g2 <- (1-R) * (cbind(Z,Xr) %*% beta1 - gamma* sigmay2) + R * Y - psi
    g <- cbind(Sy,Sz,g1,g2)
    
    return(g)
}

#####################################################################
### Estimating equations for doubly robust estimator
# Input: para = (psi, gamma, alpha, beta1, beta2, sigmay2, sigmaz2)
# Output: A n*(ds+dr+dz+5) matrix, each column stands for one estimating equation

DRmrf <-  function(para, data){
    Y <- data$Y
    R <- data$R
    Y[R==0] <- 0
    Z <- data$Z
    Xr <- data$Xr
    Xz <- data$Xz
    Xs <- data$Xs
    dr <- dim(Xr)[2]
    dz <- dim(Xz)[2]
    ds <- dim(Xs)[2]
    
    psi <- para[1]
    gamma <- para[2]
    alpha <- para[3:(ds+2)]
    beta1 <- para[(ds+3):(dr+ds+3)]
    beta2 <- para[(dr+ds+4):(dr+dz+ds+3)]
    beta12 <- beta1[1]
    sigmay2 <- para[dr+dz+ds+4]
    sigmaz2 <- para[dr+dz+ds+5]
    
    w <- INVweight(Xs,Y,R,gamma,alpha)
    beta0 <- Xz %*% beta2 - gamma * beta12 * sigmaz2
    Sr <- (w*R - 1) * cbind(Xs, Z-beta0)
    Sy <- cbind(as.numeric(Y - cbind(Z,Xr) %*% beta1) * cbind(Z,Xr),  
                (Y - cbind(Z,Xr) %*% beta1)^2 - sigmay2) / mean(R)
    Sz <- cbind(as.numeric(Z - Xz %*% beta2)*Xz, 
                (Z - Xz %*% beta2)^2 - sigmaz2) / mean(R)
    Sy[R==0,] <- 0
    Sz[R==0,] <- 0
    
    g1 <- numeric(length(Y))
    g1[R==1] <- w[R==1] * Y[R==1]
    g1 <- g1 + (1-w*R)*(cbind(Z,Xr) %*% beta1 - gamma*sigmay2) - psi
    g <- cbind(Sr,Sy,Sz,g1)
    
    return(g)
}

#####################################################################
### Estimating equations for semiparametric efficient estimator
# Input: para = (psi, gamma, alpha, beta1, beta2, sigmay2, sigmaz2)
# Output: A n*4 matrix, each column stands for one estimating equation

EFFmrf <- function(para, data){
    Y <- data$Y
    R <- data$R
    Y[R==0] <- 0
    Z <- data$Z
    Xr <- data$Xr
    Xz <- data$Xz
    Xs <- data$Xs
    X <- cbind(Xr,Xz,Xs)
    dr <- dim(Xr)[2]
    dz <- dim(Xz)[2]
    ds <- dim(Xs)[2]
    
    psi <- para[1]
    gamma <- para[2]
    alpha <- para[3:(ds+2)]
    beta1 <- para[(ds+3):(dr+ds+3)]
    beta2 <- para[(dr+ds+4):(dr+dz+ds+3)]
    beta12 <- beta1[1]
    sigmay2 <- para[dr+dz+ds+4]
    sigmaz2 <- para[dr+dz+ds+5]
    
    IntQ <- function(z, xr, xz, xs, sigmaZ, gam2sigmay2){
        meanZ = xz %*% beta2 - gamma * beta12 * sigmaz2
        dnorm(z, meanZ, sigmaZ) * plogis(c(xs%*%alpha) + gamma * (z*beta12 + c(xr%*%beta1[-1])) - 3/2*gam2sigmay2)
    }

    IntPartialQ <- function(z, xr, xz, xs, sigmaZ, gamsigmay2, gam2sigmay2){
        meanZ = c(xz %*% beta2) - gamma * beta12 * sigmaz2
        logit = plogis(c(xs%*%alpha) + gamma*(z*beta12 + c(xr%*%beta1[-1])) - 3/2*gam2sigmay2)
        parQ = logit * (1 - logit) * ((z*beta12 + c(xr%*%beta1[-1])) - 3*gamsigmay2)
        parDens = logit * beta12 * (z - (c(xz%*%beta2) - gamma*beta12*sigmaz2))
        
        dnorm(z, mean = meanZ, sd = sigmaZ) * (parQ - parDens)
    }
    
    Qzfunc <- function(z, xr, xz, xs, sigmaZ, gam2sigmay2){
        meanZ = c(xz %*% beta2) - gamma * beta12 * sigmaz2
        dnorm(z, mean = meanZ, sd = sigmaZ) * plogis(c(xs%*%alpha) + gamma*(z*beta12 + c(xr%*%beta1[-1])) - 3/2*gam2sigmay2) * z
    }
    
    w <- INVweight(Xs,Y,R,gamma,alpha)
    wR = w*R
    negwR = 1-w*R
    muy = cbind(Xz %*% beta2, Xr) %*% beta1
    sigmay = sigmay2 + (beta12^2) * sigmaz2
    gamsigmay = gamma * sigmay
    gamsigmay2 = gamma * sigmay2
    gam2sigmay2 = gamma^2 * sigmay2
    muz = cbind(Z, Xr) %*% beta1
    sigmaZ = sqrt(sigmaz2)
    beta12sigmaz2 = sigmaz2*beta12^2
    Xsalpha = Xs%*%alpha
    
    IF0 = negwR * (muy - gamsigmay) + wR*Y - psi
    
    Q = plogis(Xsalpha + gamma*muz - 3*gam2sigmay2/2)
    negQ = 1-Q
    K = muy - muz - gamma*beta12sigmaz2 + gamsigmay2*negQ
    EQ = partialEQ = EQz = numeric(length(Y))
    for (i in 1:dim(X)[1]){
      EQ[i]=integrate(IntQ, lower=-Inf, upper=Inf, xr=X[i,1:dr], xz=X[i,(dr+1):(dr+dz)], xs=X[i,(dr+dz+1):(dr+dz+ds)], sigmaZ, gam2sigmay2)$value
    }
    EK = sigmay2 * gamma * (1-EQ)
    projIF0 = negwR * (K - Q*EK/EQ)
    IF1 = IF0 - projIF0
    
    partialw = R*Y*exp(-gamma*Y - Xsalpha)
    partialIF0 = partialw * (muy - gamsigmay - Y) - negwR * sigmay
    partialQ = Q * negQ * (muz - 3*gamsigmay2)
    partialK = sigmay2 * (negQ- gamma*partialQ) - beta12sigmaz2
    
    for (i in 1:dim(X)[1]){
        partialEQ[i] = integrate(IntPartialQ, lower=-Inf, upper=Inf, xr=X[i,1:dr], xz=X[i,(dr+1):(dr+dz)], xs=X[i,(dr+dz+1):(dr+dz+ds)], sigmaZ, gamsigmay2, gam2sigmay2)$value
    }
    
    partialEK = sigmay2 * (1 - EQ - gamma*partialEQ)
    partialprojIF0 = partialw * (K - Q * EK / EQ) + negwR * (partialK - ((partialQ*EK + Q*partialEK)*EQ - Q*EK*partialEQ) / (EQ^2))
    partialIF1 = partialIF0 - partialprojIF0
    EpartialIF1 = mean(partialIF1)
    
    fr = plogis(Xsalpha + gamma*muz - gam2sigmay2 / 2)
    Sgamma = (fr - R) * (gamsigmay2 - muz)
    
    for (i in 1:dim(X)[1]){
        EQz[i] = integrate(Qzfunc, lower=-Inf, upper=Inf, xr=X[i,1:dr], xz=X[i,(dr+1):(dr+dz)], xs=X[i,(dr+dz+1):(dr+dz+ds)], sigmaZ, gam2sigmay2)$value
    }
    effSgamma = (1 - wR) * beta12 * Q * (EQz / EQ - Z)
    EIFgamma = effSgamma / mean(effSgamma^2)
    EIFpsi = IF1 + EpartialIF1 * EIFgamma

    return(cbind(EIFpsi+psi, EIFgamma+gamma))
}

EFFmrf2 <- function(para, data, EIF){
    Y <- data$Y
    R <- data$R
    Y[R==0] <- 0
    Z <- data$Z
    Xr <- data$Xr
    Xz <- data$Xz
    Xs <- data$Xs
    dr <- dim(Xr)[2]
    dz <- dim(Xz)[2]
    ds <- dim(Xs)[2]
    
    psieff <- para[1]
    gammaeff <- para[2]
    psi <- para[3]
    gamma <- para[4]
    alpha <- para[5:(ds+4)]
    beta1 <- para[(ds+5):(dr+ds+5)]
    beta2 <- para[(dr+ds+6):(dr+dz+ds+5)]
    beta12 <- beta1[1]
    sigmay2 <- para[dr+dz+ds+6]
    sigmaz2 <- para[dr+dz+ds+7]
    
    w <- INVweight(Xs,Y,R,gamma,alpha)
    beta0 <- Xz %*% beta2 - gamma * beta12 * sigmaz2
    Sr <- (w*R - 1) * cbind(Xs, Z-beta0)
    Sy <- cbind(as.numeric(Y - cbind(Z,Xr) %*% beta1) * cbind(Z,Xr),  
                (Y - cbind(Z,Xr) %*% beta1)^2 - sigmay2) / mean(R)
    Sz <- cbind(as.numeric(Z - Xz %*% beta2)*Xz, 
                (Z - Xz %*% beta2)^2 - sigmaz2) / mean(R)
    Sy[R==0,] <- 0
    Sz[R==0,] <- 0
    
    g1 <- g2 <- g3 <- numeric(length(Y))
    g1 <- EIF[,1] - psieff
    g2 <- EIF[,2] - gammaeff
    g3[R==1] <- w[R==1] * Y[R==1]
    g3 <- g3 + (1-w*R)*(cbind(Z,Xr) %*% beta1 - gamma*sigmay2) - psi
    g <- cbind(g1,g2,g3,Sr,Sy,Sz)
    
    return(g)
}

# Derivative of score equations
G1 <- function(bfun,para,data){
    G1 <- apply(bfun(para,data), 2, mean)
    return(G1)
}
G12 <- function(bfun,para,data,EIF){
    G1 <- apply(bfun(para,data,EIF), 2, mean)
    return(G1)
}

G <- function(bfun,para,data){
    G <- jacobian(func=G1,bfun=bfun,x=para,data=data)
    return(G)
}
G2 <- function(bfun,para,data,EIF){
    G <- jacobian(func=G12,bfun=bfun,x=para,data=data,EIF=EIF)
    return(G)
}

# Variance estimation
VarEst <- function(bfun,para,data){
    bG <- solve(G(bfun,para,data) + diag(10^-7,length(para)))
    bg <- bfun(para,data)
    spsz <- dim(bg)[1]
    Omega <- t(bg)%*%bg/spsz
    Sigma <- bG%*%Omega%*%t(bG)
    
    return(Sigma/spsz)
}
VarEst2 <- function(bfun,para,data,EIF){
    bG <- solve(G2(bfun,para,data,EIF) + diag(10^-7,length(para)))
    bg <- bfun(para,data,EIF)
    spsz <- dim(bg)[1]
    Omega <- t(bg)%*%bg/spsz
    Sigma <- bG%*%Omega%*%t(bG)
    
    return(Sigma/spsz)
}

# Confidence interval
ConfInt <- function(esti, ci){
    esti <- as.matrix(esti)
    dm <- dim(esti)[2]
    para <- esti[,1:(dm/2)]
    dvar <- esti[,(dm/2+1):dm]
    z <- -qnorm((1-ci)/2)
    dsd <- sqrt(dvar)
    return(list(lb=para-z*dsd, ub=para+z*dsd))
}

# Coverage probability
CoverProb <- function(esti,ci,trvlu){
    esti <- as.matrix(esti)
    dm <- dim(esti)[2]
    para <- esti[,1:(dm/2)]
    dvar <- esti[,(dm/2+1):dm]
    z <- -qnorm((1-ci)/2)
    dsd <- sqrt(dvar)
    lb <- para-z*dsd; ub <- para+z*dsd
    
    return(trvlu>=lb&trvlu<=ub)
}

# P-value based on normal approximation
Pvalue <- function(esti){
    esti <- as.matrix(esti)
    dm <- dim(esti)[2]
    para <- esti[,1:(dm/2)]
    dvar <- esti[,(dm/2+1):dm]
    dsd <- sqrt(dvar)
    return((1 - pnorm(abs(para),mean=0,sd=dsd))*2)
} 

