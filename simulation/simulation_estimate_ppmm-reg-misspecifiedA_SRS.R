########################################
# Simulation for mass imputation
# Estimates: PPMM-based method with "misspecified" model A - SRS Samples
# Author: Rebecca Andridge
# Last modified: 02/13/2026
########################################
library(tidyverse)
require(doParallel)
require(foreach)
require(doRNG)
library(survey)
library(rstanarm)

# POPNUM <- 2
# MECHNUM <- 1
# i <- 1

source("../functions/mub_reg_bayes_local.R")

# Set up parallel processing -----
cl <- makeCluster(detectCores(logical = FALSE)-1)
registerDoParallel(cl)

# CHOOSE SAMPLE SIZES:  -----
SAMPSIZES <- 1    # n_P = 5000, n_NP = 500
# SAMPSIZES <- 2    # n_P = 5000, n_NP = 10000

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
        
        # Z variables for P sample, adding column of 1s for intercept
        ZMAT_P <- as.matrix(cbind(1, samp.P[,c("Z1","A")])) # Z1 and A only, excludes Z1
        
        ## STEP 1: DRAW S|Z -------
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
          S ~ Z1 + A,
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
        # apply(draws_betas_s_model, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
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
        # check to see how often each unit's draw will come from S=1
        # table(apply(DRAWS_S, 1, sum))
        
        ## STEP 2: Draw adjusted regression coefficients from PPMM and impute Y|Z from mixture ------
        
        ### PPPM with Uniform(0,1) prior on phi -----
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
        est_ppmmregmisA <- est$qbar
        var_ppmmregmisA <- est$t
        se_ppmmregmisA <- sqrt(var_ppmmregmisA)
        tstar <- qt(0.975, est$df)
        lb_ppmmregmisA <- est_ppmmregmisA - tstar*se_ppmmregmisA
        ub_ppmmregmisA <- est_ppmmregmisA + tstar*se_ppmmregmisA
        
        ### PPPM with phi = 1 -----
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_NS.srs, zparams = 2, ndraws = NIMPS, userphi = 1)
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
        est_ppmmregmisA1 <- est$qbar
        var_ppmmregmisA1 <- est$t
        se_ppmmregmisA1 <- sqrt(var_ppmmregmisA1)
        tstar <- qt(0.975, est$df)
        lb_ppmmregmisA1 <- est_ppmmregmisA1 - tstar*se_ppmmregmisA1
        ub_ppmmregmisA1 <- est_ppmmregmisA1 + tstar*se_ppmmregmisA1

        ### PPPM with phi = 0 -----
        set.seed(SEEDS[i])
        draws <- mub_reg_bayes(stats_NP, stats_NS.srs, zparams = 2, ndraws = NIMPS, userphi = 0)
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
        est_ppmmregmisA0 <- est$qbar
        var_ppmmregmisA0 <- est$t
        se_ppmmregmisA0 <- sqrt(var_ppmmregmisA0)
        tstar <- qt(0.975, est$df)
        lb_ppmmregmisA0 <- est_ppmmregmisA0 - tstar*se_ppmmregmisA0
        ub_ppmmregmisA0 <- est_ppmmregmisA0 + tstar*se_ppmmregmisA0
        
        # Results -----
        estimates <- data.frame(cbind(est_ppmmregmisA,  se_ppmmregmisA,  lb_ppmmregmisA,  ub_ppmmregmisA,
                                      est_ppmmregmisA1, se_ppmmregmisA1, lb_ppmmregmisA1, ub_ppmmregmisA1,
                                      est_ppmmregmisA0, se_ppmmregmisA0, lb_ppmmregmisA0, ub_ppmmregmisA0))
        estimates$probtype <- "SRS"
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
    filename4 <- paste0("./simresults/ppmmregmisA_SRS_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
    save(res, file = filename4)
  }
}
stop.time <- Sys.time()
# check runtime
stop.time - start.time

# deregister parallel processing
stopCluster(cl)

# sizes = 1
# 1000 reps = Time difference of 16.674 hours (WINDOWS)
# sizes = 2
# 1000 reps =  (WINDOWS)