########################################
# Simulation for mass imputation
# Estimates: PPMM-based methods
# Author: Rebecca Andridge
# Last modified: 06/28/2025
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)

source("../functions/mub_reg_bayes_local.R")
source("../functions/mub_reg_local.R")

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
    
    set.seed(28102)
    SEEDS <- round(runif(NREPS,1,999999))
    # Run in parallel -----
    res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("ppmm")) %dopar%
      {
        # Create P and NP samples using sampled indices -----
        # P sample (not including Y)
        samp.P <- pop[Smat.P[,i],c("Z1","Z2","A")]
        n_P <- nrow(samp.P)
        # NP sample
        samp.NP <- pop[Smat.NP[,i],]
        n_NP <- nrow(samp.NP)
        
        # Inference using new prediction idea -----
        # summary statistics for NP sample (selected sample)
        stats_NP <- list(n_s = n_NP,
                         mean_YZA_s = colMeans(samp.NP),
                         var_YZA_s = var(samp.NP))
        # summary statistics for P sample (non-selected sample)
        stats_P <- list(n_ns = n_P,
                        mean_ZA_ns = colMeans(samp.P),
                        var_ZA_ns = var(samp.P))
        # Extract only Z variables of interest from P sample, adding column of 1s for intercept
        ZMAT_P <- as.matrix(cbind(1, samp.P[,c("Z1","Z2")])) # Z1 and Z2 only, excludes A
        
        ## Uniform(0,1) prior on phi
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_P, zparams = 2, ndraws = 1000)
        # predict Y in P sample (using draws of adjusted betas)
        ybarhats <- ZMAT_P %*% t(draws$betadraws)   # rows = observations, columns = draws
        residerrors <- apply(as.matrix(draws$residdraws), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        yhat <- ybarhats + residerrors
        # draws of estimated mean of Y in P sample (SRS so weights = 1)
        draws_ymean <- colMeans(yhat)
        # posterior median and 95% CrI for overall mean Y
        est_ppmmreg <- median(draws_ymean)
        # credible interval
        lb_ppmmreg <- as.numeric(quantile(draws_ymean, 0.025))
        ub_ppmmreg <- as.numeric(quantile(draws_ymean, 0.975))
        
        ## set phi = 1
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_P, zparams = 2, ndraws = 1000, userphi = 1)
        # predict Y in P sample (using draws of adjusted betas)
        ybarhats <- ZMAT_P %*% t(draws$betadraws)   # rows = observations, columns = draws
        residerrors <- apply(as.matrix(draws$residdraws), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        yhat <- ybarhats + residerrors
        # draws of estimated mean of Y in P sample (SRS so weights = 1)
        draws_ymean <- colMeans(yhat)
        # posterior median and 95% CrI for overall mean Y
        est_ppmmreg1 <- median(draws_ymean)
        # credible interval
        lb_ppmmreg1 <- as.numeric(quantile(draws_ymean, 0.025))
        ub_ppmmreg1 <- as.numeric(quantile(draws_ymean, 0.975))

        ## set phi = 0
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_P, zparams = 2, ndraws = 1000, userphi = 0)
        # predict Y in P sample (using draws of adjusted betas)
        ybarhats <- ZMAT_P %*% t(draws$betadraws)   # rows = observations, columns = draws
        residerrors <- apply(as.matrix(draws$residdraws), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        yhat <- ybarhats + residerrors
        # draws of estimated mean of Y in P sample (SRS so weights = 1)
        draws_ymean <- colMeans(yhat)
        # posterior median and 95% CrI for overall mean Y
        est_ppmmreg0 <- median(draws_ymean)
        # credible interval
        lb_ppmmreg0 <- as.numeric(quantile(draws_ymean, 0.025))
        ub_ppmmreg0 <- as.numeric(quantile(draws_ymean, 0.975))
        
        # Inference using direct application of PPMM -----
        # summary statistics for NP sample (selected sample)
        stats_NP <- list(n_YZ = n_NP,
                         mean_YZ = colMeans(samp.NP),
                         var_YZ = var(samp.NP))
        # summary statistics for P sample (non-selected sample)
        stats_P <- list(n_Z = n_P,
                        mean_Z = colMeans(samp.P),
                        var_Z = var(samp.P))

        ## Uniform prior on phi
        set.seed(SEEDS[i])
        draws <- ppmm::means_bayes(stats_NP,
                                   stats_P,
                                   phi_character = "runif(1)",
                                   ndraws = 1000)
        # draws of mean of Y in P sample
        est_ppmm <- median(draws$muY_ns)
        # 95% credible interval
        lb_ppmm <- as.numeric(quantile(draws$muY_ns, 0.025))
        ub_ppmm <- as.numeric(quantile(draws$muY_ns, 0.975))
        
        ## set phi = 1
        set.seed(SEEDS[i])
        draws <- ppmm::means_bayes(stats_NP,
                                   stats_P,
                                   phi = 1,
                                   ndraws = 1000)
        # draws of mean of Y in P sample
        est_ppmm1 <- median(draws$muY_ns)
        # 95% credible interval
        lb_ppmm1 <- as.numeric(quantile(draws$muY_ns, 0.025))
        ub_ppmm1 <- as.numeric(quantile(draws$muY_ns, 0.975))

        ## set phi = 0
        set.seed(SEEDS[i])
        draws <- ppmm::means_bayes(stats_NP,
                                   stats_P,
                                   phi = 0,
                                   ndraws = 1000)
        # draws of mean of Y in P sample
        est_ppmm0 <- median(draws$muY_ns)
        # 95% credible interval
        lb_ppmm0 <- as.numeric(quantile(draws$muY_ns, 0.025))
        ub_ppmm0 <- as.numeric(quantile(draws$muY_ns, 0.975))
        
        # Return results -----
        data.frame(cbind(est_ppmmreg, lb_ppmmreg, ub_ppmmreg,
                         est_ppmmreg1, lb_ppmmreg1, ub_ppmmreg1,
                         est_ppmmreg0, lb_ppmmreg0, ub_ppmmreg0,
                         est_ppmm, lb_ppmm, ub_ppmm,
                         est_ppmm1, lb_ppmm1, ub_ppmm1,
                         est_ppmm0, lb_ppmm0, ub_ppmm0))
      }
    # Add pop and mech numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM
    
    # Save out results -----
    filename3 <- paste0("./simresults/ppmmreg_pop",POPNUM,"mech",MECHNUM,".RData")
    save(res, file = filename3)
  }
}
  
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# Time difference of 53.4327 mins
