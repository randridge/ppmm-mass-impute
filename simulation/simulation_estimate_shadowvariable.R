########################################
# Simulation for mass imputation
# Estimates: Shadow Variable approach
# Author: Rebecca Andridge
# Last modified: 06/30/2025
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)

# POPNUM <- 3
# MECHNUM <- 3

source('../functions/BasFunctions.R')
source("../functions/BasFunctions_edit.R")

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  # Load population data -----
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Loop through selection mechanisms
  for (MECHNUM in 1:9)
  {
    print(paste0("Mech:", MECHNUM))
    # Load samples info -----
    filename2 <- paste0("./simdata/samp_pop",POPNUM,"mech",MECHNUM,".RData")
    load(filename2)
    NREPS <- ncol(Smat.NP)
    
    # Run in parallel -----
    res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("numDeriv")) %dopar%
      {
        # Create P and NP samples using sampled indices -----
        # P sample (not including Y)
        samp.P <- pop[Smat.P[,i],c("Z1","Z2","A")]
        samp.P$Y <- NA
        n_P <- nrow(samp.P)
        # NP sample
        samp.NP <- pop[Smat.NP[,i],]
        n_NP <- nrow(samp.NP)
        # stack
        both <- rbind(samp.P, samp.NP)
        
        # Shadow variable approach -----
        Y <- both$Y # outcome
        R <- !is.na(Y) # response indicator
        Z <- both$A  # shadow variable
        Xr <- Xz <- Xs <- X <- as.matrix(cbind(1,both[,c("Z1","Z2")])) # covariates plus intercept
        data1 <- list(Y=Y, R=R, Z=Z, Xs=Xs, Xr=Xr, Xz=Xz)
        # some starting values
        ymean_NP <- mean(Y, na.rm = T)
        # baseline propensity model
        fitr <- glm(R ~ X - 1, family=binomial())
        # baseline outcome estimation
        fity <- lm(Y[R==1] ~ Z[R==1] + X[R==1,] - 1)
        # ancillary outcome estimation
        fitz <- lm(Z[R==1] ~ X[R==1,] - 1)
        ### DR estimate
        # starting param values
        drstartpar <- c(ymean_NP,            # psi = outcome mean
                        1,                   # gamma = OR
                        coef(fitr),          # coef from propensity model
                        fity$coef,           # coef from Y|Z,X,R=1
                        fitz$coef,           # coef from Z|X,R=1
                        summary(fity)$sig^2, # resid variance from Y|Z,X,R=1
                        summary(fitz)$sig^2, # resid variance from Z|X,R=1
                        ymean_NP)            # psiNR = NP outcome mean
        # point estimates
        drpar <- optim(par = drstartpar, # starting values of parameter vector
                       fn = GMM,         # function to be min/max'd
                       f = DRmrf_edit, data=data1, # additional arguments for fn 
                       method = "BFGS", hessian = FALSE)$par
        # variance estimates
        drvar <- VarEst(DRmrf_edit, drpar, data1)
        # overall mean Y
        est_shadowdr <- as.numeric(drpar[1])
        # SE of mean Y
        se_shadowdr <- sqrt(drvar[1,1])
        # 95% CI
        z <- qnorm(.975)
        lb_shadowdr <- est_shadowdr - z*se_shadowdr
        ub_shadowdr <- est_shadowdr + z*se_shadowdr
        # back-calculate mean in probability sample (easy since SRS)
        est_shadowdrP <- ((n_NP+n_P)*est_shadowdr - n_NP*ymean_NP)/n_P
        # note that drpar[15] = ymean_NP
        # SE of mean Y in prob sample
        se_shadowdrP <- sqrt((n_NP+n_P)^2*drvar[1,1] + n_NP^2*drvar[15,15] - 2*(n_NP+n_P)*n_NP*drvar[1,15])/n_P
        # 95% CI
        lb_shadowdrP <- est_shadowdrP - z*se_shadowdrP
        ub_shadowdrP <- est_shadowdrP + z*se_shadowdrP
        # estimated OR(Y)
        or_shadowdr <- drpar[2]
        orse_shadowdr <- sqrt(drvar[2,2])

        # Return results -----
        data.frame(cbind(est_shadowdr, se_shadowdr, lb_shadowdr, ub_shadowdr,
                         est_shadowdrP, se_shadowdrP, lb_shadowdrP, ub_shadowdrP,
                         or_shadowdr, orse_shadowdr))
      }
   # Add pop and mech numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM

    # Save out results -----
    filename3 <- paste0("./simresults/shadowvar_pop",POPNUM,"mech",MECHNUM,".RData")
    save(res, file = filename3)
  }
}
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# Time difference of 2.267163 hours
