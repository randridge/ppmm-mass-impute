#####################################################################
### Estimating equations for doubly robust estimator
# Input: para = (psi, gamma, alpha, beta1, beta2, sigmay2, sigmaz2)
# Output: A n*(ds+dr+dz+5) matrix, each column stands for one estimating equation

DRmrf_edit <-  function(para, data){
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
  
  # add estimating equation for mean for R=1
  psiR <- para[dr+dz+ds+6]
  g2 <- Y - psiR
  g2[R==0] <- 0
  
  g <- cbind(Sr,Sy,Sz,g1,g2)
  
  return(g)
}