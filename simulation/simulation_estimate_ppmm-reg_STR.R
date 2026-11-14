########################################
# Simulation for mass imputation
# Estimates: PPMM-based method - Stratified Samples
# Author: Rebecca Andridge
# Last modified: 02/10/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(survey)
library(rstanarm)

# POPNUM <- 3
# MECHNUM <- 9
# i <- 1

source("../functions/mub_reg_bayes_local.R")

# Set up parallel processing -----
cl <- makeCluster(detectCores(logical = FALSE)-1)
registerDoParallel(cl)

# CHOOSE SAMPLE SIZES:  -----
# SAMPSIZES <- 1    # n_P = 5000, n_NP = 500
SAMPSIZES <- 2    # n_P = 5000, n_NP = 10000

NREPS <- 1000

# Loop through populations 
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
    
    set.seed(28102)
    SEEDS <- round(runif(NREPS,1,999999))
    
# START <- Sys.time()
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
        
        # Create sample using sampled indices
        samp.P <- pop[Smat.P.STRAT[,i],]
        samp.P$WT <- wts$wt_strat
        samp.P$N <- ifelse(samp.P$STRATUM == 1, N1, N2)
        n_P <- nrow(samp.P)
        des <- survey::svydesign(~1, weights = ~WT, strata = ~STRATUM, fpc = ~N, data = samp.P)
        mean_P <- coef(survey::svymean(~Z1+Z2+A, des))
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
        stats_NS.str <- list(n_ns = n_P,
                             mean_ZA_ns = xbar_ns,
                             var_ZA_ns = s2_ns)
        # stratum sample sizes for later calculations
        n1 <- sum(samp.P$STRATUM == 1)
        n2 <- nrow(samp.P) - n1
        
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
          iter = 1500,) # 500 draws after warm-up
          # cores = parallel::detectCores())
        
        # posterior draws of logistic regression coefficients
        draws_betas_s_model <- as.matrix(bayesmodel_s)
        # apply(draws_betas_s_model, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
        # thin the chain for use in imputation to only 100 draws
        draws_betas_s_model <- draws_betas_s_model[seq(5,500,5),]
        
        # check summary
        # summary(bayesmodel_s)
        
        # check trace plots
        # plot(bayesmodel_s, plotfun = "trace")
        
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
        # check to see how often each unit's draw will come from S=1
        # table(apply(DRAWS_S, 1, sum))
        
        # STEP 2: Draw adjusted regression coefficients from PPMM and impute Y|Z from mixture ------
        ## PPPM with Uniform(0,1) prior on phi -----
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_NS.str, zparams = 2, ndraws = NIMPS)
        # Draws of Y|Z,S=1 in P sample
        DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Draws of Y|Z,S=0 in P sample
        DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Combine to get draws of Y|Z (mixture distribution)
        DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
        # Design-based estimates of mean of y in each of the NIMPS imputed datasets
        # using hand-coded calculations b/c much faster than survey::svymean()
        theta <- apply(DRAW_y, 2, function(Y) sum(Y*samp.P$WT)/sum(samp.P$WT))
        var_within <- apply(DRAW_y, 2, function(Y) {(N1/N)^2 * var(Y[samp.P$STRATUM == 1])/n1 * (1 - n1/N1) + (N2/N)^2 * var(Y[samp.P$STRATUM == 2])/n2 * (1 - n2/N2)})
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
        draws <- mub_reg_bayes(stats_NP, stats_NS.str, zparams = 2, ndraws = NIMPS, userphi = 0)
        # Draws of Y|Z,S=1 in P sample
        DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Draws of Y|Z,S=0 in P sample
        DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Combine to get draws of Y|Z (mixture distribution)
        DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
        # Design-based estimates of mean of y in each of the NIMPS imputed datasets
        theta <- apply(DRAW_y, 2, function(Y) sum(Y*samp.P$WT)/sum(samp.P$WT))
        var_within <- apply(DRAW_y, 2, function(Y) {(N1/N)^2 * var(Y[samp.P$STRATUM == 1])/n1 * (1 - n1/N1) + (N2/N)^2 * var(Y[samp.P$STRATUM == 2])/n2 * (1 - n2/N2)})
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
        draws <- mub_reg_bayes(stats_NP, stats_NS.str, zparams = 2, ndraws = NIMPS, userphi = 1)
        # Draws of Y|Z,S=1 in P sample
        DRAW_y_s1 <- ZMAT_P %*% t(draws$betadraws_s) + apply(as.matrix(draws$residdraws_s), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Draws of Y|Z,S=0 in P sample
        DRAW_y_s0 <- ZMAT_P %*% t(draws$betadraws_ns) + apply(as.matrix(draws$residdraws_ns), 1, function(x) rnorm(n_P, mean = 0, sd = sqrt(x)))
        # Combine to get draws of Y|Z (mixture distribution)
        DRAW_y <- DRAW_y_s1*DRAWS_S + DRAW_y_s0*(1-DRAWS_S)
        # Design-based estimates of mean of y in each of the NIMPS imputed datasets
        theta <- apply(DRAW_y, 2, function(Y) sum(Y*samp.P$WT)/sum(samp.P$WT))
        var_within <- apply(DRAW_y, 2, function(Y) {(N1/N)^2 * var(Y[samp.P$STRATUM == 1])/n1 * (1 - n1/N1) + (N2/N)^2 * var(Y[samp.P$STRATUM == 2])/n2 * (1 - n2/N2)})
        # combine with Rubin's Rules
        est <- mice::pool.scalar(Q = theta, U = var_within, n = n_P)
        est_ppmmreg1 <- est$qbar
        var_ppmmreg1 <- est$t
        se_ppmmreg1 <- sqrt(var_ppmmreg1)
        tstar <- qt(0.975, est$df)
        lb_ppmmreg1 <- est_ppmmreg1 - tstar*se_ppmmreg1
        ub_ppmmreg1 <- est_ppmmreg1 + tstar*se_ppmmreg1
        
        # Results -----
        estimates <- data.frame(cbind(est_ppmmreg,  se_ppmmreg,  lb_ppmmreg,  ub_ppmmreg,
                                      est_ppmmreg1, se_ppmmreg1, lb_ppmmreg1, ub_ppmmreg1,
                                      est_ppmmreg0, se_ppmmreg0, lb_ppmmreg0, ub_ppmmreg0))
        estimates$probtype <- "STR"
        estimates$rep <- i
        estimates$num_s_equals_1 <- sum(DRAWS_prob_s_raw/(1-DRAWS_prob_s_raw) > 1)
        return(estimates)
      }
# STOP <- Sys.time()
    # Add pop and selection mechanism and sample sizes numbers -----
    res$pop <- POPNUM
    res$mech <- MECHNUM
    res$sampsize <- SAMPSIZES
    
    # Save out results -----
    filename4 <- paste0("./simresults/ppmmreg_STR_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
    save(res, file = filename4)
  }
}
stop.time <- Sys.time()
# check runtime
stop.time - start.time

# deregister parallel processing
stopCluster(cl)

# sizes = 1
# all 27 * 1000 reps = Time difference of 14.33106 hours (Windows)
# sizes = 2
# all 27 * 1000 reps = Time difference of 1.569873 days (Windows)
