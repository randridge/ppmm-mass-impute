####################################
# Generate data and samples for simulation
# Author: Rebecca R Andridge
# Last Modified: 01/17/2026
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

# Find intercept in selection model for NP sample -----
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
# FIRST SET: Sampling fraction = 500/100000 = 0.005
parms1 <- as.data.frame(t(mapply(function(a,b,c,d,e,f,g,h) findg0(0.005,a,b,c,d,e,f,g,h), 
                                 a=sim.parms$rho_y1, 
                                 b=sim.parms$rho_y2, 
                                 c=sim.parms$rho_1a, 
                                 d=sim.parms$rho_ya.z, 
                                 e=sim.parms$gamma_y, 
                                 f=sim.parms$gamma_1, 
                                 g=sim.parms$gamma_2, 
                                 h=sim.parms$gamma_a)))
# SECOND SET: Sampling fraction = 10000/100000 = 0.1
parms2 <- as.data.frame(t(mapply(function(a,b,c,d,e,f,g,h) findg0(0.1,a,b,c,d,e,f,g,h), 
                                 a=sim.parms$rho_y1, 
                                 b=sim.parms$rho_y2, 
                                 c=sim.parms$rho_1a, 
                                 d=sim.parms$rho_ya.z, 
                                 e=sim.parms$gamma_y, 
                                 f=sim.parms$gamma_1, 
                                 g=sim.parms$gamma_2, 
                                 h=sim.parms$gamma_a)))
# stack
parms <- rbind(parms1, parms2)
names(parms) <- c("sfrac", 
                  "rho_y1", "rho_y2", "rho_1a", "rho_ya.z", "rho_ya", 
                  "gamma_0", "gamma_y", "gamma_1", "gamma_2", "gamma_a")
rm(parms1, parms2)

# Add odds ratios for selection model
parms$or_y <- exp(parms$gamma_y)
parms$or_1 <- exp(parms$gamma_1)
parms$or_2 <- exp(parms$gamma_2)
parms$or_a <- exp(parms$gamma_a)

# Number the populations and NP selection mechanisms and relative sample sizes
# 3 populations
sub <- parms %>% 
  dplyr::select(rho_y1, rho_y2, rho_1a, rho_ya.z) %>% 
  distinct() %>% 
  mutate(pop = 1:n())
# 9 NP selection mechanisms
sub2 <- parms %>% 
  dplyr::select(gamma_y, gamma_1, gamma_2, gamma_a) %>% 
  distinct() %>%
  mutate(mech = 1:n())
# 2 sets of sample sizes
sub3 <- parms %>%
  dplyr::select(sfrac) %>% 
  distinct() %>%
  mutate(sampsize = 1:n())
# combine
parms <- full_join(sub, parms)
parms <- full_join(sub2, parms)
parms <- full_join(sub3, parms)
rm(sub, sub2, sub3, sim.parms, findg0)

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
  # stratum labels based on Z1
  pop$STRATUM <- ifelse(pop$Z1 < 1, 1, 2)
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  save(pop, file=filename)
}

for (i in 1:3)
{
  gendata(i, 531217)
}
rm(i, gendata)

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

rm(truth, POPNUM, filename)

# Draw probability samples -----
drawPsamps <- function(POPNUM, NREPS, SEED)
{
  # Load population data
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  
  # Sample size
  N <- nrow(pop)
  n <- 5000
  
  # Indices for units in each stratum, for easier sampling
  indices_stratum1 <- c(1:N)[pop$STRATUM == 1]
  indices_stratum2 <- c(1:N)[pop$STRATUM == 2]
  
  # Number of units in each stratum
  N1 <- sum(pop$STRATUM == 1)
  N2 <- sum(pop$STRATUM == 2)
  
  # Set seed
  set.seed(SEED)
  
  # Draw samples NREPS times
  Smat.P.SRS <- matrix(nrow = n, ncol = NREPS)
  Smat.P.STRAT <- matrix(nrow = n, ncol = NREPS)
  for (i in 1:NREPS)
  {
    # Draw SRSWOR sample
    Smat.P.SRS[,i] <- sample(1:N, n, replace = FALSE)
    
    # Draw stratified sample with equal allocation, total sample size = n_P
    Smat.P.STRAT[,i] <- c(sample(indices_stratum1, n/2, replace = FALSE),  # stratum 1
                          sample(indices_stratum2, n/2, replace = FALSE))  # stratum 2
  }
  
  # Sampling weights
  wts <- list(wt_srs = N/n,
              wt_strat = c(rep(N1/(n/2), n/2), rep(N2/(n/2), n/2)))
  
  # Save matrices with indices of selected units and weights info
  filename2 <- paste0("./simdata/Psamp_pop",POPNUM,".RData")
  
  save(Smat.P.SRS, Smat.P.STRAT, wts, file=filename2)
}

for (pop in 1:3)
{
  print(paste("Population:", pop))
  drawPsamps(pop, 1000, 2025999)
}
rm(pop)

# Draw nonprobability samples -----
drawNPsamps <- function(POPNUM, MECHNUM, SAMPSIZES, NREPS, SEED)
{
  # Load simulation parameters and subset to scenario being considered
  load("./simdata/simParms.RData")
  pop.parms <- subset(parms, pop==POPNUM & mech==MECHNUM & sampsize==SAMPSIZES)
  
  # Load population data and drop the stratum indicator
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  pop$STRATUM <- NULL
  
  # Population size
  N <- nrow(pop)
  
  # Sample size
  if (SAMPSIZES == 1)
  {
    n <- 500
  } else if (SAMPSIZES == 2)
  {
    n <- 10000
  }
  
  # Set seed
  set.seed(SEED)
  
  # Draw samples NREPS times
  Smat.NP <- matrix(nrow = n, ncol = NREPS)
  for (i in 1:NREPS)
  {
    # Calculate probabilities of selection
    linpred <- drop(cbind(rep(1,N),as.matrix(pop)) %*% t(pop.parms[c("gamma_0", "gamma_y", "gamma_1", "gamma_2", "gamma_a")]))
    prob <- plogis(linpred)
    # Draw sample (without replacement) with specified unequal probs to get sample of desired size
    pik <- sampling::inclusionprobabilities(prob, n)
    sind.NP <- sampling::UPsystematic(pik, eps = 1e-07)
    Smat.NP[,i] <- c(1:N)[sind.NP == 1]
  }
  
  # Save matrix with indices of selected units
  filename2 <- paste0("./simdata/NPsamp_pop",POPNUM,"mech",MECHNUM,"relsizes",SAMPSIZES,".RData")
  save(Smat.NP, file=filename2)
}

for (pop in 1:3)
{
  print(paste("Population:", pop))
  for (mech in 1:9)
  {
    print(paste("Mechanism:", mech))
    for (sampsizes in 1:2)
    {
      print(paste("Relative Sample Sizes:", sampsizes))
      drawNPsamps(pop, mech, sampsizes, 1000, 2025999)
    }
  }
}

# check selection probabilities -----
load("./simdata/simParms.RData")

get_probrange <- function(POPNUM, MECHNUM, SAMPSIZES)
{
  pop.parms <- subset(parms, pop==POPNUM & mech==MECHNUM & sampsize==SAMPSIZES)
  
  # Load population data and drop the stratum indicator
  filename <- paste0("./simdata/pop",POPNUM,".RData")
  load(filename)
  pop$STRATUM <- NULL
  
  # Population size
  N <- nrow(pop)
  
  # Sample size
  if (SAMPSIZES == 1)
  {
    n <- 500
  } else if (SAMPSIZES == 2)
  {
    n <- 10000
  }
  
  # Calculate probabilities of selection
  linpred <- drop(cbind(rep(1,N),as.matrix(pop)) %*% t(pop.parms[c("gamma_0", "gamma_y", "gamma_1", "gamma_2", "gamma_a")]))
  prob <- plogis(linpred)
  pik <- sampling::inclusionprobabilities(prob, n)
  # Summary stats for probs
  return(list(probs = pik,
              sumstats = c(pop = POPNUM, 
                           mech = MECHNUM, 
                           sampsize = SAMPSIZES, 
                           min_pik = min(pik),
                           max_pik = max(pik),
                           mean_pik = mean(pik),
                           median_pik = median(pik),
                           p75_pik = as.numeric(quantile(pik, 0.75)),
                           p90_pik = as.numeric(quantile(pik, 0.90)),
                           p95_pik = as.numeric(quantile(pik, 0.95)))))
}

results <- vector()
probs <- vector()
for (POPNUM in 1:3)
{
  for (MECHNUM in 1:9)
  {
    for (SAMPSIZES in 1:2)
    {
      res <- get_probrange(POPNUM, MECHNUM, SAMPSIZES)
      results <- rbind(results, res$sumstats)
      temp <- tibble(pop = POPNUM, mech = MECHNUM, sampsize = SAMPSIZES, prob = res$probs)
      probs <- rbind(probs, temp)
    }
  }
}
rm(temp)
results <- as_tibble(results) %>%
  full_join(parms)
results <- results %>%
  mutate(mech_f = factor(mech, levels = 1:9,
                         labels = c("OR(Y) = 1, OR(A) = 1",
                                    "OR(Y) = 1.1, OR(A) = 1",
                                    "OR(Y) = 2, OR(A) = 1",
                                    "OR(Y) = 1, OR(A) = 1.1",
                                    "OR(Y) = 1.1, OR(A) = 1.1",
                                    "OR(Y) = 2, OR(A) = 1.1",
                                    "OR(Y) = 1, OR(A) = 2",
                                    "OR(Y) = 1.1, OR(A) = 2",
                                    "OR(Y) = 2, OR(A) = 2")))
theme_set(theme_bw())

pd <- position_dodge(0.5)
# range of probs
# results %>%
#   ggplot(aes(y = mech_f, x = max_pik)) +
#   geom_errorbar(aes(xmin = min_pik, xmax = max_pik), width = 0.2) +
#   facet_grid(sampsize ~ rho_ya.z, labeller = label_both) + 
#   theme(legend.position = "top")
results %>%
  ggplot(aes(y = mech_f, x = max_pik, col = factor(rho_ya.z))) +
  geom_errorbar(aes(xmin = min_pik, xmax = max_pik), width = 0.5, position = pd) +
  facet_grid(sampsize ~ ., labeller = label_both) + 
  theme(legend.position = "top")
ggsave("./simulation_prob_ranges.pdf", width = 8, height = 8)
write_csv(results, "simulation_prob_ranges.csv")

# mean, median of probs
results %>%
  filter(sampsize == 1) %>%
  ggplot(aes(y = mech_f, x = mean_pik)) +
  geom_point(aes(x = median_pik), color = "red", shape = 1) +
  geom_point(aes(x = p75_pik), color = "orange", shape = 3) +
  geom_point(aes(x = p95_pik), color = "green", shape = 4) +
  facet_grid(sampsize ~ rho_ya.z, labeller = label_both) + 
  theme(legend.position = "top")
ggsave("./simulation_prob_mean_median.pdf", width = 8, height = 6)


probs <- probs %>%
  left_join(parms)
probs <- probs %>%
  mutate(mech_f = factor(mech, levels = 1:9,
                         labels = c("OR(Y) = 1, OR(A) = 1",
                                    "OR(Y) = 1.1, OR(A) = 1",
                                    "OR(Y) = 2, OR(A) = 1",
                                    "OR(Y) = 1, OR(A) = 1.1",
                                    "OR(Y) = 1.1, OR(A) = 1.1",
                                    "OR(Y) = 2, OR(A) = 1.1",
                                    "OR(Y) = 1, OR(A) = 2",
                                    "OR(Y) = 1.1, OR(A) = 2",
                                    "OR(Y) = 2, OR(A) = 2")))
probs %>%
  ggplot(aes(y = mech_f, x = prob, col = factor(rho_ya.z), fill = factor(rho_ya.z))) +
  geom_violin(position = pd) +
  facet_grid(sampsize ~ ., labeller = label_both) + 
  theme(legend.position = "top")
ggsave("./simulation_prob_violin.pdf", width = 8, height = 8)
