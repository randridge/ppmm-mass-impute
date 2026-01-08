####################################
# Generate data and samples for simulation
# Author: Rebecca R Andridge
# Last Modified: 06/27/2025
####################################

require(tidyverse)

# Set model parms -----
sim.parms <- expand.grid(rho_y1=0.4, 
                         rho_y2=0.4, 
                         rho_1a=0.4, 
                         rho_ya.z=c(0.2,0.5,0.8), 
                         gamma_y=c(0, log(1.1), log(2)), 
                         gamma_1=log(1.1), 
                         gamma_2=log(1.1), 
                         gamma_a=c(0, log(1.1), log(2)))

# Find intercept in selection model -----
# sfrac=0.05; rho_y1=0.6; rho_y2=0.6; rho_1a=0.6; rho_ya.z=0.8; gamma_y=log(2); gamma_1=log(2); gamma_2=log(2); gamma_a=log(2)
findg0 <- function(sfrac, rho_y1, rho_y2, rho_1a, rho_ya.z, gamma_y, gamma_1, gamma_2, gamma_a)
{
  set.seed(531217)
  N <- 100000
  ####
  # Fixed for the simulation
  ####
  mu_y <- 10                 # E[Y] = 10
  mu_1 <- mu_2 <- mu_a <- 0  # E[Z1] = E[Z2] = E[A] = 0
  s_yy <- 4                  # Var(Y) = 4
  s_11 <- s_22 <- s_aa <- 1  # Var(Z1) = Var(Z2) = Var(A) = 1
  s_12 <- 0                  # Cov(Z1,Z2) = 0
  s_2a <- 0                  # Cov(Z2,A) = 0
  ####
  # Calculate covariances based on input correlations
  ####
  s_y1 <- rho_y1*sqrt(s_yy*s_11)   # Cov(Y,Z1)
  s_y2 <- rho_y2*sqrt(s_yy*s_22)   # Cov(Y,Z2)
  s_1a <- rho_1a*sqrt(s_11*s_aa)   # Cov(Z1,A)
  denom <- (s_yy-(s_y1^2+s_y2^2))*(s_aa-(s_1a^2+s_2a^2))
  if (denom > 0) {
    s_ya <- rho_ya.z*sqrt(denom) + s_y1*s_1a + s_y2*s_2a  # Cov(Y,A)
    ####
    # Calculate unconditional correlation Corr(Y,A)
    ####
    rho_ya <- s_ya/sqrt(s_yy*s_aa)
    ####
    # Create covariance matrix
    ####
    SIGMA <- matrix(c(s_yy, s_y1, s_y2, s_ya,
                      s_y1, s_11, s_12, s_1a,
                      s_y2, s_12, s_22, s_2a,
                      s_ya, s_1a, s_2a, s_aa), nrow=4, byrow=TRUE)
    ####
    # Generate dataset
    ####
    YZA <- mnormt::rmnorm(N, mean=c(mu_y, mu_1, mu_2, mu_a), varcov=SIGMA)

    ####
    # Find value of intercept for specified selection fraction
    ####
    f <- function(gamma_0) abs(sfrac - mean(1 - 1/(1+exp(drop(gamma_0 + YZA %*% c(gamma_y, gamma_1, gamma_2, gamma_a))))))
    gamma_0 <- optimize(f, c(-20, 0))$minimum
  } else {
    rho_ya <- gamma_0  <- NA
  }
  return(c(sfrac, 
           rho_y1, rho_y2, rho_1a, rho_ya.z, rho_ya, 
           gamma_0, gamma_y, gamma_1, gamma_2, gamma_a))
}

# Calculate the intercept and estimate the selection bias for each set of simulation parameters
# Sampling fraction = 500/100000 = 0.005
parms <- as.data.frame(t(mapply(function(a,b,c,d,e,f,g,h) findg0(0.005,a,b,c,d,e,f,g,h), 
                                a=sim.parms$rho_y1, 
                                b=sim.parms$rho_y2, 
                                c=sim.parms$rho_1a, 
                                d=sim.parms$rho_ya.z, 
                                e=sim.parms$gamma_y, 
                                f=sim.parms$gamma_1, 
                                g=sim.parms$gamma_2, 
                                h=sim.parms$gamma_a)))
names(parms) <- c("sfrac", 
                  "rho_y1", "rho_y2", "rho_1a", "rho_ya.z", "rho_ya", 
                  "gamma_0", "gamma_y", "gamma_1", "gamma_2", "gamma_a")

# Add odds ratios for selection model
parms$or_y <- exp(parms$gamma_y)
parms$or_1 <- exp(parms$gamma_1)
parms$or_2 <- exp(parms$gamma_2)
parms$or_a <- exp(parms$gamma_a)

# Number the populations and selection mechanisms
sub <- parms %>% 
  dplyr::select(rho_y1, rho_y2, rho_1a, rho_ya.z) %>% 
  distinct() %>% 
  mutate(pop = 1:n())
sub2 <- parms %>% 
  dplyr::select(gamma_y, gamma_1, gamma_2, gamma_a) %>% 
  distinct() %>%
  mutate(mech = 1:n())
parms <- full_join(sub, parms)
parms <- full_join(sub2, parms)
rm(sub, sub2,sim.parms, findg0)

# Generate population datasets -----
gendata <- function(POPNUM, SEED)
{
  # Fixed for the simulation
  N <- 100000                # Population size
  mu_y <- 10                 # E[Y] = 10
  mu_1 <- mu_2 <- mu_a <- 0  # E[Z1] = E[Z2] = E[A] = 0
  s_yy <- 4                  # Var(Y) = 4
  s_11 <- s_22 <- s_aa <- 1  # Var(Z1) = Var(Z2) = Var(A) = 1
  s_12 <- 0                  # Cov(Z1,Z2) = 0
  s_2a <- 0                  # Cov(Z2,A) = 0
  
  # Pull off additional parameters for the chosen population/mechanism
  pop.parms <- subset(parms, pop==POPNUM & mech==1) # doesn't matter which selection mech
  
  # Calculate covariances based on input correlations
  s_y1 <- pop.parms$rho_y1[1]*sqrt(s_yy*s_11)   # Cov(Y,Z1)
  s_y2 <- pop.parms$rho_y2[1]*sqrt(s_yy*s_22)   # Cov(Y,Z2)
  s_1a <- pop.parms$rho_1a[1]*sqrt(s_11*s_aa)   # Cov(Z1,A)
  s_ya <- pop.parms$rho_ya[1]*sqrt(s_yy*s_aa)   # Cov(Y,A)
  
  # Set seed 
  set.seed(SEED)
  
  # Generate population data
  YZA <- mnormt::rmnorm(N, mean=c(mu_y, mu_1, mu_2, mu_a), varcov=matrix(c(s_yy, s_y1, s_y2, s_ya,
                                                                           s_y1, s_11, s_12, s_1a,
                                                                           s_y2, s_12, s_22, s_2a,
                                                                           s_ya, s_1a, s_2a, s_aa), nrow=4, byrow=TRUE))
  pop <- as.data.frame(YZA)
  names(pop) <- c("Y","Z1","Z2","A")
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  save(pop, file=filename)
}

for (i in 1:3)
{
  gendata(i, 531217)
}
rm(i, gendata)

# Draw samples -----
sampdata <- function(POPNUM, MECHNUM, NREPS, SEED)
{
  # Pull off selection model parameters for the chosen population/mechanism
  pop.parms <- subset(parms, pop==POPNUM & mech==MECHNUM)
  
  # Load population data
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Sample sizes
  N <- nrow(pop)            # pop size
  n_NP <- N*pop.parms$sfrac # NP sample
  n_P <- 10*n_NP            # P sample
  
  # Set seed
  set.seed(SEED)
  
  # Draw samples NREPS times
  Smat.P <- matrix(nrow = n_P, ncol = NREPS)
  Smat.NP <- matrix(nrow = n_NP, ncol = NREPS)
  for (i in 1:NREPS)
  {
    # Sample selection: PROBABILITY SAMPLE
    # Draw SRSWOR sample
    Smat.P[,i] <- sample(1:N, n_P, replace = FALSE)
    
    # Sample selection: NONPROBABILITY SAMPLE
    # Calculate probabilities of selection
    linpred <- drop(cbind(rep(1,N),as.matrix(pop)) %*% t(pop.parms[c("gamma_0", "gamma_y", "gamma_1", "gamma_2", "gamma_a")]))
    prob <- plogis(linpred)
    # Draw sample (without replacement) with specified unequal probs to get sample of desired size
    pik <- sampling::inclusionprobabilities(prob, n_NP)
    sind.NP <- sampling::UPsystematic(pik, eps = 1e-07)
    Smat.NP[,i] <- c(1:N)[sind.NP == 1]
  }
  # save matrices with indices of selected units
  filename2 <- paste0("./simdata/samp_pop",POPNUM,"mech",MECHNUM,".RData")
  
  save(Smat.P, Smat.NP, file=filename2)
}

for (pop in 1:3)
{
  for (mech in 1:9)
  {
    sampdata(pop, mech, 1000, 2025999)
  }
}

# Get true pop values -----
truth <- matrix(nrow = 3, ncol = 5)

for (POPNUM in 1:3)
{
  print(paste0("Pop:", POPNUM))
  
  # Load population data -----
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  truth[POPNUM,1] <- POPNUM
  
  # True mean of Y -----
  truth[POPNUM,2] <- mean(pop$Y)
  
  # Regression parameter estimates from Y|Z -----
  truth[POPNUM,3:5] <- as.numeric(coef(lm(Y ~ Z1 + Z2, data = pop)))
}

truth <- as.data.frame(truth)
names(truth) <- c("pop", "mu", "B0", "B1", "B2")

# add to parms dataset
parms <- left_join(parms, truth)

# save parms
save(parms, file="./simdata/simParms.RData")


##### ADDITIONAL BELOW ##### -----
library(ggplot2)
theme_set(theme_bw())

# Difference between selected sample regression coefficients and pop regression coefficients

parms$b0_s_bias <- 100*(parms$b0_s - parms$B0)/parms$B0
parms$b1_s_bias <- 100*(parms$b1_s - parms$B1)/parms$B1
parms$b2_s_bias <- 100*(parms$b2_s - parms$B2)/parms$B2

ggplot(parms, aes(x=factor(or_a), y=b0_s_bias, color=factor(rho_ya.z))) + 
  geom_point(position=position_dodge(0.5)) + 
  facet_grid(~or_y, labeller = label_both)

ggplot(parms, aes(x=factor(or_a), y=b1_s_bias, color=factor(rho_ya.z))) + 
  geom_point(position=position_dodge(0.5)) + 
  facet_grid(~or_y, labeller = label_both)

ggplot(parms, aes(x=factor(or_a), y=b2_s_bias, color=factor(rho_ya.z))) + 
  geom_point(position=position_dodge(0.5)) + 
  facet_grid(~or_y, labeller = label_both)

