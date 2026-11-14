########################################
# Simulation for shadow variable approach
# Estimates: Chen DR/MI/IPW approach
# Author: Rebecca Andridge
# Last modified: 02/05/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(nonprobsvy)

# POPNUM <- 3
# MECHNUM <- 9

source('../functions/BasFunctions.R')
source("../functions/BasFunctions_edit.R")

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# CHOOSE SAMPLE SIZES:  -----
# SAMPSIZES <- 1    # n_P = 5000, n_NP = 500
SAMPSIZES <- 2    # n_P = 5000, n_NP = 10000

NREPS <- 1000

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  # Load population data -----
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Loop through NP selection mechanisms
  for (MECHNUM in 1:9)
  {
    print(paste0("Selection:", MECHNUM))
    
    # Load matrices with lists of sampled units - NP samples
    filename2 <- paste0("./simdata/NPsamp_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
    load(filename2)
    # Load matrices with lists of sampled units - P samples
    filename3 <- paste0("./simdata/Psamp_pop",POPNUM,".RData")
    load(filename3)
# START <- Sys.time()    
    # Run in parallel -----
    res_try <- foreach(i=1:NREPS, .inorder=FALSE, .packages = c("numDeriv"),
                   .errorhandling = "pass") %dopar%
      {
        # P sample: SRS -----
        # Create sample using sampled indices
        samp.srs <- pop[Smat.P.SRS[,i],c("Z1","Z2","A")]
        samp.srs$Y <- NA
        n_P <- nrow(samp.srs)
        
        # NP sample -----
        # Create sample using sampled indices
        samp.NP <- pop[Smat.NP[,i],c("Y","Z1","Z2","A")]
        n_NP <- nrow(samp.NP)
        
        # stack
        both <- rbind(samp.srs, samp.NP)
        
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
        # overall mean Y (NP + P samples together)
        est_shadowdr <- as.numeric(drpar[1])
        # back-calculate mean in probability sample (easy since SRS)
        est_shadowdrP <- ((n_NP+n_P)*est_shadowdr - n_NP*ymean_NP)/n_P
        # note that drpar[15] = ymean_NP
        # SE of mean Y in prob sample
        se_shadowdrP <- sqrt((n_NP+n_P)^2*drvar[1,1] + n_NP^2*drvar[15,15] - 2*(n_NP+n_P)*n_NP*drvar[1,15])/n_P
        # 95% CI
        z <- qnorm(.975)
        lb_shadowdrP <- est_shadowdrP - z*se_shadowdrP
        ub_shadowdrP <- est_shadowdrP + z*se_shadowdrP
        # estimated OR(Y)
        or_shadowdr <- drpar[2]
        orse_shadowdr <- as.numeric(sqrt(drvar[2,2]))
        
        # Return results -----
        estimates <- data.frame(cbind(est_shadowdrP, se_shadowdrP, lb_shadowdrP, ub_shadowdrP,
                                      or_shadowdr, orse_shadowdr))
        
        estimates$probtype <- "SRS"
        estimates$rep <- i
        
        estimates
      }
# STOP <- Sys.time()
    # only keep replicates that worked
    res_list <- purrr::keep(res_try, is.data.frame)
    res <- list_rbind(res_list)

    # Add pop and selection mechanism and sample sizes numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM
    res$sampsize <- SAMPSIZES
      
    # Save out results -----
    filename4 <- paste0("./simresults/shadowvar_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
    save(res, file = filename4)
  }
}
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# MAC: 1000 reps for 1 POP / 1 MECH / smaller NP sample: Time difference of 3.525326 mins
# MAC: 1000 reps for 2 POP / 9 MECH / smaller NP sample: Time difference of 1.264949 hours
# MAC 100 reps / 1 POP / 9 MECHS / larger NP sample: Time difference of 13.92107 mins
# WINDOWS 100 reps / 2 POPs / 9 MECH / larger NP sample: Time difference of 2.207694 hours
# WINDOWS 1000 reps / 2 POPs / 9 MECH / larger NP sample: Time difference of 22.23044 hours
# WINDOWS 1000 reps / 1 POPs / 9 MECH / larger NP sample: Time difference of 10.84893 hours