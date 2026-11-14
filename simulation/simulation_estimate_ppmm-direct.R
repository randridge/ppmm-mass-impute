########################################
# Simulation for mass imputation
# Estimates: Direct application of PPMM
# Author: Rebecca Andridge
# Last modified: 01/18/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(survey)
library(ppmm)

# POPNUM <- 1
# MECHNUM <- 1
# SAMPSIZES <- 2

# Set up parallel processing -----
cl <- makeCluster(detectCores(logical = FALSE)-1)
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
  N <- N1 + N2
  
  # Loop through relative sample sizes -----
  for (SAMPSIZES in 1:2)
  {
    print(paste0("RelSize:", SAMPSIZES))

    # Loop through selection mechanisms
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
      
      set.seed(28102)
      SEEDS <- round(runif(NREPS,1,999999))

      # Run in parallel -----
      res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE, .packages = c("survey", "ppmm")) %dopar%
        {
          # NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],c("Y","Z1","Z2","A")]
          n_NP <- nrow(samp.NP)
          stats_NP <- list(n_YZ = n_NP,
                           mean_YZ = colMeans(samp.NP),
                           var_YZ = var(samp.NP))
          
          # P sample by SRS -----
          # Create sample using sampled indices
          samp.srs <- pop[Smat.P.SRS[,i],]
          samp.srs$WT <- wts$wt_srs
          samp.srs$N <- N1 + N2
          n_P.srs <- nrow(samp.srs)
          des.srs <- survey::svydesign(~1, weights = ~WT, fpc = ~N, data = samp.srs)
          mean_P.srs <- coef(survey::svymean(~Z1+Z2+A, des.srs))
          var_P.srs <- as.matrix(survey::svyvar(~Z1+Z2+A, des.srs))[1:3,1:3]
          # back out the selected estimates to get nonselected estimates
          xbar <- mean_P.srs
          s2 <- var_P.srs
          xbar_s <- stats_NP$mean_YZ[-1]
          s2_s <- stats_NP$var_YZ[-1,-1]
          sfrac <- n_NP/N
          n_s <- n_NP
          n_ns <- N - n_s
          # back out the means
          xbar_ns <- (xbar - sfrac*xbar_s) / (1-sfrac)
          # back out the variance matrix
          s2_ns <- ((n_s + n_ns - 1) * s2 - (n_s - 1) * s2_s - n_s * xbar_s %*% t(xbar_s) - n_ns * xbar_ns %*% t(xbar_ns) + (n_s + n_ns) * xbar %*% t(xbar))/(n_ns-1)
          # backed out NS sample stats
          stats_NS.srs <- list(n_Z = n_ns,
                               mean_Z = xbar_ns,
                               var_Z = s2_ns)
          
          # P sample: Stratified -----
          # Create sample using sampled indices
          samp.str <- pop[Smat.P.STRAT[,i],]
          samp.str$WT <- wts$wt_strat
          samp.str$N <- ifelse(samp.str$STRATUM == 1, N1, N2)
          n_P.str <- nrow(samp.str)
          des.str <- survey::svydesign(~1, weights = ~WT, strata = ~STRATUM, fpc = ~N, data = samp.str)
          mean_P.str <- coef(survey::svymean(~Z1+Z2+A, des.str))
          var_P.str <- as.matrix(survey::svyvar(~Z1+Z2+A, des.str))[1:3,1:3]
          # back out the selected estimates to get nonselected estimates
          xbar <- mean_P.str
          s2 <- var_P.str
          xbar_s <- stats_NP$mean_YZ[-1]
          s2_s <- stats_NP$var_YZ[-1,-1]
          sfrac <- n_NP/N
          n_s <- n_NP
          n_ns <- N - n_s
          # back out the means
          xbar_ns <- (xbar - sfrac*xbar_s) / (1-sfrac)
          # back out the variance matrix
          s2_ns <- ((n_s + n_ns - 1) * s2 - (n_s - 1) * s2_s - n_s * xbar_s %*% t(xbar_s) - n_ns * xbar_ns %*% t(xbar_ns) + (n_s + n_ns) * xbar %*% t(xbar))/(n_ns-1)
          # backed out NS sample stats
          stats_NS.str <- list(n_Z = n_ns,
                               mean_Z = xbar_ns,
                               var_Z = s2_ns)
          
          # Inference using direct application of PPMM -----

          ##### Uniform prior on phi
          ### SRS
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.srs,
                                     phi_character = "runif(1)",
                                     ndraws = 1000)
          # draws of mean of Y in population (selected + nonselected)
          est_ppmm_srs <- median(draws$muY)
          se_ppmm_srs <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm_srs <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm_srs <- as.numeric(quantile(draws$muY, 0.975))
          ### Stratified
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.str,
                                     phi_character = "runif(1)",
                                     ndraws = 1000)
          # draws of mean of Y in P sample
          est_ppmm_str <- median(draws$muY)
          se_ppmm_str <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm_str <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm_str <- as.numeric(quantile(draws$muY, 0.975))
          
          ##### set phi = 1
          ### SRS
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.srs,
                                     phi = 1,
                                     ndraws = 1000)
          # draws of mean of Y in population (selected + nonselected)
          est_ppmm1_srs <- median(draws$muY)
          se_ppmm1_srs <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm1_srs <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm1_srs <- as.numeric(quantile(draws$muY, 0.975))
          ### STRATIFIED
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.str,
                                     phi = 1,
                                     ndraws = 1000)
          # draws of mean of Y in P sample
          est_ppmm1_str <- median(draws$muY)
          se_ppmm1_str <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm1_str <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm1_str <- as.numeric(quantile(draws$muY, 0.975))
          
          ##### set phi = 0
          ### SRS
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.srs,
                                     phi = 0,
                                     ndraws = 1000)
          # draws of mean of Y in population (selected + nonselected)
          est_ppmm0_srs <- median(draws$muY)
          se_ppmm0_srs <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm0_srs <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm0_srs <- as.numeric(quantile(draws$muY, 0.975))
          ### STRATIFIED
          set.seed(SEEDS[i])
          draws <- ppmm::means_bayes(stats_NP,
                                     stats_NS.str,
                                     phi = 0,
                                     ndraws = 1000)
          # draws of mean of Y in P sample
          est_ppmm0_str <- median(draws$muY)
          se_ppmm0_str <- sd(draws$muY)
          # 95% credible interval
          lb_ppmm0_str <- as.numeric(quantile(draws$muY, 0.025))
          ub_ppmm0_str <- as.numeric(quantile(draws$muY, 0.975))
          
          # Return results -----
          estimates <- data.frame(rbind(c(est_ppmm_srs,  se_ppmm_srs,  lb_ppmm_srs,  ub_ppmm_srs,
                                          est_ppmm1_srs, se_ppmm1_srs, lb_ppmm1_srs, ub_ppmm1_srs,
                                          est_ppmm0_srs, se_ppmm0_srs, lb_ppmm0_srs, ub_ppmm0_srs),
                                        c(est_ppmm_str,  se_ppmm_str,  lb_ppmm_str,  ub_ppmm_str,
                                          est_ppmm1_str, se_ppmm1_str, lb_ppmm1_str, ub_ppmm1_str,
                                          est_ppmm0_str, se_ppmm0_str, lb_ppmm0_str, ub_ppmm0_str)))
          names(estimates) <- c("est_ppmm",  "se_ppmm",  "lb_ppmm",  "ub_ppmm",
                                "est_ppmm1", "se_ppmm1", "lb_ppmm1", "ub_ppmm1",
                                "est_ppmm0", "se_ppmm0", "lb_ppmm0", "ub_ppmm0")
          estimates$probtype <- c("SRS", "STR")
          estimates$rep <- i
          
          return(estimates)
        }
     
      # Add pop and selection mechanism and sample sizes numbers -----
      res$pop <- POPNUM
      res$mech <- MECHNUM
      res$sampsize <- SAMPSIZES
      
      # Save out results -----
      filename4 <- paste0("./simresults/ppmmdirect_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      save(res, file = filename4)
    }
  }
}

stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# 1000 replicates, 1 population, 2 sizes, 9 mechanisms (Mac): Time difference of 10.90959 mins

# res <- as_tibble(res)
# res %>%
#   group_by(probtype) %>%
#   summarize(n = n(),
#             nreps = n(),
#             mean_est_ppmm = mean(est_ppmm),
#             empse_est_ppmm = sd(est_ppmm),
#             avgse_est_ppmm = mean(se_ppmm),
#             mean_est_ppmm0 = mean(est_ppmm0),
#             empse_est_ppmm0 = sd(est_ppmm0),
#             avgse_est_ppmm0 = mean(se_ppmm0))
