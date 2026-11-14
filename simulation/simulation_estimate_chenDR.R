########################################
# Simulation for mass imputation
# Estimates: Chen DR/MI/IPW approach
# Author: Rebecca Andridge
# Last modified: 01/17/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(nonprobsvy)

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
  
  # Stratum population sizes
  N1 <- sum(pop$STRATUM == 1)
  N2 <- sum(pop$STRATUM == 2)
  
  # Loop through relative sample sizes -----
  for (SAMPSIZES in 1:2)
  {
    print(paste0("RelSize:", SAMPSIZES))
    
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
      
      NREPS <- ncol(Smat.NP)
      
      # Run in parallel -----
      res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("survey", "nonprobsvy")) %dopar%
        {
          # P sample: SRS -----
          # Create sample using sampled indices
          samp.srs <- pop[Smat.P.SRS[,i],]
          samp.srs$WT <- wts$wt_srs
          samp.srs$N <- N1 + N2
          des.srs <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.srs)

          # P sample: Stratified -----
          # Create sample using sampled indices
          samp.str <- pop[Smat.P.STRAT[,i],]
          samp.str$WT <- wts$wt_strat
          samp.str$N <- ifelse(samp.str$STRATUM == 1, N1, N2)
          des.str <- survey::svydesign(~1, weights = ~WT, strata = ~STRATUM, fpc = ~N, data = samp.str)
          
          # NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],]
          n_NP <- nrow(samp.NP)
          
          # Nonprobsvy methods -----
          # Doubly robust
          chen_dr_srs <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                             outcome = Y ~ Z1 + Z2 + A, 
                                             data = samp.NP, 
                                             svydesign = des.srs, 
                                             method_outcome = "glm", 
                                             family_outcome = "gaussian")
          chen_dr_str <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                             outcome = Y ~ Z1 + Z2 + A, 
                                             data = samp.NP, 
                                             svydesign = des.str, 
                                             method_outcome = "glm", 
                                             family_outcome = "gaussian")
          
          # Mass imputation
          chen_mi_srs <- nonprobsvy::nonprob(outcome = Y ~ Z1 + Z2 + A, 
                                             data = samp.NP, 
                                             svydesign = des.srs, 
                                             method_outcome = "glm", 
                                             family_outcome = "gaussian")
          chen_mi_str <- nonprobsvy::nonprob(outcome = Y ~ Z1 + Z2 + A, 
                                             data = samp.NP, 
                                             svydesign = des.str, 
                                             method_outcome = "glm", 
                                             family_outcome = "gaussian")
          
          # IPW
          chen_ipw_srs <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                              target = ~ Y, 
                                              data = samp.NP, 
                                              svydesign = des.srs)
          chen_ipw_str <- nonprobsvy::nonprob(selection = ~ Z1 + Z2 + A, 
                                              target = ~ Y, 
                                              data = samp.NP, 
                                              svydesign = des.str)

          # estimates from each method/prob sample type
          # mean
          est_chendr_srs  <- as.numeric(chen_dr_srs$output[1])
          est_chenmi_srs  <- as.numeric(chen_mi_srs$output[1])
          est_chenipw_srs <- as.numeric(chen_ipw_srs$output[1])
          est_chendr_str  <- as.numeric(chen_dr_str$output[1])
          est_chenmi_str  <- as.numeric(chen_mi_str$output[1])
          est_chenipw_str <- as.numeric(chen_ipw_str$output[1])
          # standard errors
          se_chendr_srs  <- as.numeric(chen_dr_srs$output[2])
          se_chenmi_srs  <- as.numeric(chen_mi_srs$output[2])
          se_chenipw_srs <- as.numeric(chen_ipw_srs$output[2])
          se_chendr_str  <- as.numeric(chen_dr_str$output[2])
          se_chenmi_str  <- as.numeric(chen_mi_str$output[2])
          se_chenipw_str <- as.numeric(chen_ipw_str$output[2])
          # 95% CI
          lb_chendr_srs  <- as.numeric(chen_dr_srs$confidence_interval[1])
          ub_chendr_srs  <- as.numeric(chen_dr_srs$confidence_interval[2])
          lb_chenmi_srs  <- as.numeric(chen_mi_srs$confidence_interval[1])
          ub_chenmi_srs  <- as.numeric(chen_mi_srs$confidence_interval[2])
          lb_chenipw_srs <- as.numeric(chen_ipw_srs$confidence_interval[1])
          ub_chenipw_srs <- as.numeric(chen_ipw_srs$confidence_interval[2])
          lb_chendr_str  <- as.numeric(chen_dr_str$confidence_interval[1])
          ub_chendr_str  <- as.numeric(chen_dr_str$confidence_interval[2])
          lb_chenmi_str  <- as.numeric(chen_mi_str$confidence_interval[1])
          ub_chenmi_str  <- as.numeric(chen_mi_str$confidence_interval[2])
          lb_chenipw_str <- as.numeric(chen_ipw_str$confidence_interval[1])
          ub_chenipw_str <- as.numeric(chen_ipw_str$confidence_interval[2])
          
          
          # Return results -----
          estimates <- data.frame(rbind(c(est_chendr_srs,  se_chendr_srs,  lb_chendr_srs,  ub_chendr_srs,
                                          est_chenmi_srs,  se_chenmi_srs,  lb_chenmi_srs,  ub_chenmi_srs,
                                          est_chenipw_srs, se_chenipw_srs, lb_chenipw_srs, ub_chenipw_srs),
                                        c(est_chendr_str,  se_chendr_str,  lb_chendr_str,  ub_chendr_str,
                                          est_chenmi_str,  se_chenmi_str,  lb_chenmi_str,  ub_chenmi_str,
                                          est_chenipw_str, se_chenipw_str, lb_chenipw_str, ub_chenipw_str)))
          
          names(estimates) <- c("est_chendr",  "se_chendr",  "lb_chendr",  "ub_chendr",
                                "est_chenmi",  "se_chenmi",  "lb_chenmi",  "ub_chenmi",
                                "est_chenipw", "se_chenipw", "lb_chenipw", "ub_chenipw")
          estimates$probtype <- c("SRS", "STR")
          estimates$rep <- i
          
          estimates
        }
      
      # Add pop and selection mechanism and sample sizes numbers -----
      res$pop <- POPNUM
      res$mech <- MECHNUM
      res$sampsize <- SAMPSIZES
      
      # Save out results -----
      filename4 <- paste0("./simresults/chen_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      save(res, file = filename4)
    }
  }
}
stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# Time difference of 12.09256 mins
