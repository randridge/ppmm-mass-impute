########################################
# Simulation for mass imputation
# Estimate the proxy strength for various PPMM-Reg models
# Author: Rebecca Andridge
# Last modified: 01/17/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(survey)

source("../functions/mub_reg_cor.R")

sumstats <- vector()

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  
  # Load population data
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  N <- nrow(pop)
  
  # Loop through relative sample sizes -----
  for (SAMPSIZES in 2)
  {
    print(paste0("RelSize:", SAMPSIZES))
    
    # Loop through NP selection mechanisms
    for (MECHNUM in 1:9)
    {
      print(paste0("Selection:", MECHNUM))
      
      # Load matrices with lists of sampled units
      filename2 <- paste0("./simdata/NPsamp_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      load(filename2)
      # Load matrices with lists of sampled units - P samples
      filename3 <- paste0("./simdata/Psamp_pop",POPNUM,".RData")
      load(filename3)
      
      NREPS <- ncol(Smat.NP)
      
      # Run in parallel -----
      res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("survey")) %dopar%
        {
          # ORIGINAL PPMM-REG SPECIFICATION -----
          ## NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],c("Y","Z1","Z2","A")]
          n_NP <- nrow(samp.NP)
          stats_NP <- list(n_s = n_NP,
                           mean_YZA_s = colMeans(samp.NP),
                           var_YZA_s = var(samp.NP))
          ## P sample: SRS -----
          # Create sample using sampled indices
          samp.P <- pop[Smat.P.SRS[,i],]
          samp.P$WT <- wts$wt_srs
          samp.P$N <- N
          n_P <- nrow(samp.P)
          des <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.P)
          mean_P <- coef(survey::svymean(~Z1+Z2+A, design = des))
          var_P <- as.matrix(survey::svyvar(~Z1+Z2+A, des))[1:3,1:3]
          # back out the selected estimates to get nonselected estimates
          xbar <- mean_P
          s2 <- var_P
          xbar_s <- stats_NP$mean_YZA_s[-1]
          s2_s <- stats_NP$var_YZA_s[-1,-1]
          sfrac <- n_NP/N
          n_s <- n_NP
          n_ns <- N - n_s
          # back out the means
          xbar_ns <- (xbar - sfrac*xbar_s) / (1-sfrac)
          # back out the variance matrix
          s2_ns <- ((n_s + n_ns - 1) * s2 - (n_s - 1) * s2_s - n_s * xbar_s %*% t(xbar_s) - n_ns * xbar_ns %*% t(xbar_ns) + (n_s + n_ns) * xbar %*% t(xbar))/(n_ns-1)
          # backed out NS sample stats
          stats_NS.srs <- list(n_ns = N - n_NP,
                               mean_ZA_ns = xbar_ns,
                               var_ZA_ns = s2_ns)
          ## MLE of proxy strength -----
          rhohat <- mub_reg_cor(stats_NP, stats_NS.srs, nZvars = 2)
          
          
          # MISSPECIFICATION 1: Use Z2 to make proxy, use Z1 and A to impute -----
          ## NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],c("Y","Z1","A","Z2")]
          n_NP <- nrow(samp.NP)
          stats_NP <- list(n_s = n_NP,
                           mean_YZA_s = colMeans(samp.NP),
                           var_YZA_s = var(samp.NP))
          
          ## P sample: SRS -----
          # Create sample using sampled indices
          samp.P <- pop[Smat.P.SRS[,i],]
          samp.P$WT <- wts$wt_srs
          samp.P$N <- N
          n_P <- nrow(samp.P)
          des <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.P)
          mean_P <- coef(survey::svymean(~Z1+A+Z2, design = des))
          var_P <- as.matrix(survey::svyvar(~Z1+A+Z2, des))[1:3,1:3]
          # back out the selected estimates to get nonselected estimates
          xbar <- mean_P
          s2 <- var_P
          xbar_s <- stats_NP$mean_YZA_s[-1]
          s2_s <- stats_NP$var_YZA_s[-1,-1]
          sfrac <- n_NP/N
          n_s <- n_NP
          n_ns <- N - n_s
          # back out the means
          xbar_ns <- (xbar - sfrac*xbar_s) / (1-sfrac)
          # back out the variance matrix
          s2_ns <- ((n_s + n_ns - 1) * s2 - (n_s - 1) * s2_s - n_s * xbar_s %*% t(xbar_s) - n_ns * xbar_ns %*% t(xbar_ns) + (n_s + n_ns) * xbar %*% t(xbar))/(n_ns-1)
          # backed out NS sample stats
          stats_NS.srs <- list(n_ns = N - n_NP,
                               mean_ZA_ns = xbar_ns,
                               var_ZA_ns = s2_ns)
          ## MLE of proxy strength -----
          rhohat_misA <- mub_reg_cor(stats_NP, stats_NS.srs, nZvars = 2)
          
          # MISSPECIFICATION 2: Use Z2 to make proxy, use Z1 to impute (ignore A) -----
          ## NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],c("Y","Z1","Z2")]
          n_NP <- nrow(samp.NP)
          stats_NP <- list(n_s = n_NP,
                           mean_YZA_s = colMeans(samp.NP),
                           var_YZA_s = var(samp.NP))
          ## P sample: SRS -----
          # Create sample using sampled indices
          samp.P <- pop[Smat.P.SRS[,i],]
          samp.P$WT <- wts$wt_srs
          samp.P$N <- N
          n_P <- nrow(samp.P)
          des <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.P)
          mean_P <- coef(survey::svymean(~Z1+Z2, design = des))
          var_P <- as.matrix(survey::svyvar(~Z1+Z2, des))[1:2,1:2]
          # back out the selected estimates to get nonselected estimates
          xbar <- mean_P
          s2 <- var_P
          xbar_s <- stats_NP$mean_YZA_s[-1]
          s2_s <- stats_NP$var_YZA_s[-1,-1]
          sfrac <- n_NP/N
          n_s <- n_NP
          n_ns <- N - n_s
          # back out the means
          xbar_ns <- (xbar - sfrac*xbar_s) / (1-sfrac)
          # back out the variance matrix
          s2_ns <- ((n_s + n_ns - 1) * s2 - (n_s - 1) * s2_s - n_s * xbar_s %*% t(xbar_s) - n_ns * xbar_ns %*% t(xbar_ns) + (n_s + n_ns) * xbar %*% t(xbar))/(n_ns-1)
          # backed out NS sample stats
          stats_NS.srs <- list(n_ns = N - n_NP,
                               mean_ZA_ns = xbar_ns,
                               var_ZA_ns = s2_ns)
          ## MLE of proxy strength -----
          rhohat_misB <- mub_reg_cor(stats_NP, stats_NS.srs, nZvars = 1)
          
          # Return estimates -----
          c(rhohat, rhohat_misA, rhohat_misB)
        }
      
      # Add pop and selection mechanism and sample sizes numbers -----
      sumstats <- rbind(sumstats, c(pop = POPNUM, mech = MECHNUM, sampsize = SAMPSIZES, 
                                    rhohat_orig = mean(res[,1]),
                                    rhohat_misA = mean(res[,2]),
                                    rhohat_misB = mean(res[,3])))
    }
  }
}
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time

rhohat <- as_tibble(sumstats)
save(rhohat, file="./simresults/rhohat_ppmmreg.RData")

# Time difference of 2.652086 mins (MAC)