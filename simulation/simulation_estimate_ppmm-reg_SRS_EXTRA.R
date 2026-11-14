########################################
# Simulation for mass imputation
# Estimates for Extra Scenarios (w/interactions)
# Author: Rebecca Andridge
# Last modified: 02/16/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(survey)
library(rstanarm)

source("../functions/mub_reg_bayes_local.R")

# Set up parallel processing -----
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# POPNUM <- 1
# MECHNUM <- 10
# SAMPSIZES <- 1

# Loop through populations -----
start.time <- Sys.time()
for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  
  # Load population data
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # pop size
  N <- nrow(pop)
  
  # Loop through relative sample sizes -----
  for (SAMPSIZES in 1)
  {
    print(paste0("RelSize:", SAMPSIZES))
    
    # Loop through NP selection mechanisms
    for (MECHNUM in 10:11)
    {
      print(paste0("Selection:", MECHNUM))
      
      # Load matrices with lists of sampled units
      filename2 <- paste0("./simdata/NPsamp_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      load(filename2)
      # Load matrices with lists of sampled units - P samples
      filename3 <- paste0("./simdata/Psamp_pop",POPNUM,".RData")
      load(filename3)
      
      NREPS <- ncol(Smat.NP)
      
      set.seed(28102)
      SEEDS <- round(runif(NREPS,1,999999))
      
      # Run in parallel -----
      res <- foreach(i=1:NREPS, .combine=rbind, .inorder=FALSE,
                     .packages = c("survey", "rstanarm", "mice")) %dopar%
        {
          NIMPS <- 100 # hard-coded in the bayesian logistic regression
          
          # NP sample -----
          # Create sample using sampled indices
          samp.NP <- pop[Smat.NP[,i],c("Y","Z1","Z2","A")]
          n_NP <- nrow(samp.NP)
          stats_NP <- list(n_s = n_NP,
                           mean_YZA_s = colMeans(samp.NP),
                           var_YZA_s = var(samp.NP))
          
          # (Naive Estimate) Estimates from nonprobability sample -----
          est_NP <- mean(samp.NP$Y)
          se_NP <- sd(samp.NP$Y)/sqrt(n_NP)
          ci <- confint(lm(Y~1, samp.NP))
          lb_NP <- ci[1]
          ub_NP <- ci[2]
          
          # P sample: SRS -----
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
          
          # Z variables for P sample, adding column of 1s for intercept
          ZMAT_P <- as.matrix(cbind(1, samp.P[,c("Z1","Z2")])) # Z1 and Z2 only, excludes A
          
          # STEP 1: DRAW S|Z -------
          # Bayesian weighted logistic regression for S|Z
          # assign weight = 1 for NP sample units
          samp.NP$WT <- 1
          # indicator for which sample (1 = NP, 0 = P)
          samp.NP$S <- 1
          samp.P$S <- 0
          # stack NP and P datasets
          stacked <- rbind(samp.NP[,c("Z1","Z2","A","WT","S")], samp.P[,c("Z1","Z2","A","WT","S")])
          # Bayes weighted logistic regression using rstan
          bayesmodel_s <- stan_glm(
            S ~ Z1 + Z2,
            data = stacked,
            family = binomial(link = "logit"),
            weights = WT,
            refresh = 0, # don't print progress
            seed = SEEDS[i],
            chains = 1,
            warmup = 1000,
            iter = 1500) # 500 draws after warm-up  
          # cores = parallel::detectCores())
          
          # posterior draws of logistic regression coefficients
          draws_betas_s_model <- as.matrix(bayesmodel_s)
          # thin the chain for use in imputation to only 100 draws
          draws_betas_s_model <- draws_betas_s_model[seq(5,500,5),]
          
          # Transform into draws of the propensities
          # Get draws of propensities P(S=1|Z)
          xb <- ZMAT_P %*% t(draws_betas_s_model)  # rows = observations, columns = draws
          DRAWS_prob_s_raw <- 1 - 1/(1 + exp(xb))
          # bias correction for large sampling fractions for NP as in https://pubmed.ncbi.nlm.nih.gov/34219260/
          DRAWS_prob_s <- DRAWS_prob_s_raw/(1-DRAWS_prob_s_raw) 
          # if prob > 0.5, this will be > 1 --> set to 1
          DRAWS_prob_s[DRAWS_prob_s>1] <- 1
          
          # Draws of S|Z (for probability sample only)
          S_vec <- rbinom(n = n_P * NIMPS, size = 1, prob = as.vector(DRAWS_prob_s))
          DRAWS_S <- matrix(S_vec, nrow = n_P, ncol = NIMPS)
          
          # STEP 2: Draw adjusted regression coefficients from PPMM and impute Y|Z from mixture ------
          ## PPPM with Uniform(0,1) prior on phi -----
          set.seed(SEEDS[i])
          draws <- mub_reg_bayes(stats_NP, stats_NS.srs, zparams = 2, ndraws = NIMPS)
          # Draws of Y|Z,S=1 in P sample
          DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Draws of Y|Z,S=0 in P sample
          DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Combine to get draws of Y|Z (mixture distribution)
          DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
          # Design-based estimates of mean of y in each of the NIMPS imputed datasets
          # using hand-coded calculations b/c much faster than survey::svymean()
          theta <- colMeans(DRAW_y)
          var_within <- apply(DRAW_y, 2, function(Y) var(Y)/n_P) * (1 - n_P/N)
          # combine with Rubin's Rules
          est <- mice::pool.scalar(Q = theta, U = var_within, n = n_P)
          est_ppmmreg <- est$qbar
          var_ppmmreg <- est$t
          se_ppmmreg <- sqrt(var_ppmmreg)
          tstar <- qt(0.975, est$df)
          lb_ppmmreg <- est_ppmmreg - tstar*se_ppmmreg
          ub_ppmmreg <- est_ppmmreg + tstar*se_ppmmreg
          
          ## PPPM with phi = 0 -----
          set.seed(SEEDS[i])
          draws <- mub_reg_bayes(stats_NP, stats_NS.srs, zparams = 2, ndraws = NIMPS, userphi = 0)
          # Draws of Y|Z,S=1 in P sample
          DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Draws of Y|Z,S=0 in P sample
          DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Combine to get draws of Y|Z (mixture distribution)
          DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
          # Design-based estimates of mean of y in each of the NIMPS imputed datasets
          theta <- colMeans(DRAW_y)
          var_within <- apply(DRAW_y, 2, function(Y) var(Y)/n_P) * (1 - n_P/N)
          # combine with Rubin's Rules
          est <- mice::pool.scalar(Q = theta, U = var_within, n = n_P)
          est_ppmmreg0 <- est$qbar
          var_ppmmreg0 <- est$t
          se_ppmmreg0 <- sqrt(var_ppmmreg0)
          tstar <- qt(0.975, est$df)
          lb_ppmmreg0 <- est_ppmmreg0 - tstar*se_ppmmreg0
          ub_ppmmreg0 <- est_ppmmreg0 + tstar*se_ppmmreg0
          
          ## PPPM with phi = 1 -----
          set.seed(SEEDS[i])
          draws <- mub_reg_bayes(stats_NP, stats_NS.srs, zparams = 2, ndraws = NIMPS, userphi = 1)
          # Draws of Y|Z,S=1 in P sample
          DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Draws of Y|Z,S=0 in P sample
          DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
          # Combine to get draws of Y|Z (mixture distribution)
          DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
          # Design-based estimates of mean of y in each of the NIMPS imputed datasets
          theta <- colMeans(DRAW_y)
          var_within <- apply(DRAW_y, 2, function(Y) var(Y)/n_P) * (1 - n_P/N)
          # combine with Rubin's Rules
          est <- mice::pool.scalar(Q = theta, U = var_within, n = n_P)
          est_ppmmreg1 <- est$qbar
          var_ppmmreg1 <- est$t
          se_ppmmreg1 <- sqrt(var_ppmmreg1)
          tstar <- qt(0.975, est$df)
          lb_ppmmreg1 <- est_ppmmreg1 - tstar*se_ppmmreg1
          ub_ppmmreg1 <- est_ppmmreg1 + tstar*se_ppmmreg1
          
          
          # Return results -----
          estimates <- data.frame(cbind(est_NP, se_NP, lb_NP, ub_NP,
                                        est_ppmmreg,  se_ppmmreg,  lb_ppmmreg,  ub_ppmmreg,
                                        est_ppmmreg1, se_ppmmreg1, lb_ppmmreg1, ub_ppmmreg1,
                                        est_ppmmreg0, se_ppmmreg0, lb_ppmmreg0, ub_ppmmreg0))
          estimates$probtype <- "SRS"
          estimates$rep <- i
          estimates$num_s_equals_1 <- sum(DRAWS_prob_s_raw/(1-DRAWS_prob_s_raw) > 1)
          
          estimates
        }
      
      # Add pop and selection mechanism and sample sizes numbers -----
      res$pop <- POPNUM
      res$mech <- MECHNUM
      res$sampsize <- SAMPSIZES
      
      # Save out results -----
      filename4 <- paste0("./simresults/ppmmreg_SRS_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
      save(res, file = filename4)
    }
  }
}

stop.time <- Sys.time()
# Deregister parallel processing
stopCluster(cl)
# check runtime
stop.time - start.time
# MAC, 100 reps: Time difference of 7.079397 mins
# WINDOWS, 1000 reps: Time difference of 3.440208 hours