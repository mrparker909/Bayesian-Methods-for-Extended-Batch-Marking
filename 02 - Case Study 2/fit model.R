library(nimble)

######################################################
#' DEFINITIONS:
#' 
#' OBSERVED DATA:
#' To  the number of capture–recapture occasions.
#' G   the number of batch-marked release groups; G ≤ T.
#' Sg  the number of individuals sampled at sampling occasion g; g = 1, 2, . . . , G.
#' Rg  the number of individuals marked and released at sampling occasion g from batch group g; g = 1, 2, . . . , G. We condition on the {Rg} when we form the likelihood for marked individuals.
#' mgt the number of individuals from batch group g recaptured at recapture occasion t; g = 1, 2, . . . , G, t = g + 1, . . . , T.
#' ut  the number of individuals captured at sampling occasion t that were not marked; t = 1, . . . , T.
#' lt the number of individuals lost at sampling occasion t; t = 2, . . . , T. ℓt = ut − Rt.
#'
#' LATENT VARIABLES:
#' Mgt the number of individuals present at occasion t from marked group g; g = 1, 2, . . . , G, t = g + 1, . . . , T.
#' dgt the number of individuals from marked group g that die at sampling occasion t. 
#' Ut  the total unmarked individuals present at occasion t
#' At  the total unmarked individuals which survive from t to t+1
#' Ct  the total unmarked individuals recruited from t to t+1
#' 
#' PARAMETERS:
#' phi    survival
#' p      capture probability
#' lambda initial unmarked population size
#' eta    recruitment rate
#'
######################################################

################################################################
# SETTING UP DATA
# Sand Lance data
counts = readxl::read_xls(path = "02 - Case Study 2/Data/SandLanceData.xls", sheet = 2)
groupData = readxl::read_xls(path = "02 - Case Study 2/Data/SandLanceData.xls", sheet=1)
cov_data = readxl::read_xls(path = "02 - Case Study 2/Data/SandLanceData.xls", sheet=3)
countsMat = as.matrix(counts)
deltaT = groupData$DeltaDays[-1]/7 # number of weeks between sampling occassions: leads to weekly survival estimates and weekly recruitments
To = ncol(countsMat)
G  = nrow(countsMat)
Rg = groupData$`Number Marked (Rg)`[-length(groupData$`Number Marked (Rg)`)]
ut = groupData$`Number Unmarked (ut)` 
lt = rep(0, times=To)
Sg = groupData$`Number Sampled (Sg)`

mgt = array(0, dim=c(G,To))
for(g in 1:G) {
  # number of individuals from group g recaptured at recapture occasion t
  mgt[g,] = c(0,countsMat[g,-1])
}

temperature = as.numeric(groupData$Temperature)
temperature = (temperature - mean(temperature, na.rm=T)) / sd(temperature, na.rm=T)
temperature.inits = rep(NA, times=length(temperature))
temperature.inits[which(is.na(temperature))] <- mean(temperature, na.rm = T)

salinity = as.numeric(groupData$Salinity)
salinity = (salinity - mean(salinity, na.rm=T)) / sd(salinity, na.rm=T)
salinity.inits = rep(NA, times=length(salinity))
salinity.inits[which(is.na(salinity))] <- mean(salinity, na.rm = T)

cov_dat_g = as.numeric(cov_data$`Batch (g)`)
cov_dat_t = as.numeric(cov_data$`Sampling Occassion (t)`)

i_weight = as.numeric(cov_data$Mass)
i_weight.inits = rep(NA, times=length(i_weight))
i_weight.inits[which(is.na(i_weight))] <- mean(i_weight, na.rm = T)

weightCounts = table(cov_dat_g, cov_dat_t)
weightCounter = 0 * weightCounts + 1
weightSize = c(To,To,max(weightCounts))
weightArray = array(NA, dim=weightSize)
for(i in 1:length(i_weight)) {
  loc = c(cov_dat_g[i], cov_dat_t[i])
  weightArray[loc[1], loc[2], weightCounter[loc[1], loc[2]]] = i_weight[i]
  weightCounter[loc[1], loc[2]] = weightCounter[loc[1], loc[2]] + 1
}
weight.inits = weightArray * NA
weight.inits[which(is.na(weightArray), arr.ind = T)] <- mean(weightArray, na.rm = T)

i_length = as.numeric(cov_data$Length)
i_length.inits = rep(NA, times=length(i_length))
i_length.inits[which(is.na(i_length))] <- mean(i_length, na.rm = T)

lengthCounts = table(cov_dat_g, cov_dat_t)
lengthCounter = 0 * lengthCounts + 1
lengthSize = c(To,To,max(lengthCounts))
lengthArray = array(NA, dim=lengthSize)
for(i in 1:length(i_length)) {
  loc = c(cov_dat_g[i], cov_dat_t[i])
  lengthArray[loc[1], loc[2], lengthCounter[loc[1], loc[2]]] = i_length[i]
  lengthCounter[loc[1], loc[2]] = lengthCounter[loc[1], loc[2]] + 1
}

# standardize lengths:
lengthArray <- (lengthArray - mean(lengthArray, na.rm = T)) / sd(lengthArray, na.rm = T)

# we need to reorder NA's to the end of the list
meanLengths = matrix(0, nrow=To, ncol=To)
for(g in 1:To) {
  for(t in 1:To) {
    lengthArray[g,t,] = c(lengthArray[g,t,!is.na(lengthArray[g,t,])], lengthArray[g,t,is.na(lengthArray[g,t,])])
    meanLengths[g,t] = mean(lengthArray[g,t,], na.rm = T)
  }
}
meanLengths[which(is.na(meanLengths))] <- mean(lengthArray, na.rm = T)

length.inits = lengthArray * NA
length.inits[which(is.na(lengthArray), arr.ind = T)] <- mean(lengthArray, na.rm = T)

# percent missing data
1 - length(lengthArray[which(!is.na(lengthArray))]) / (7*sum(diag(lengthCounts)))

# max individuals at g,t with lengths measured
maxI = lengthCounts

# DONE SETTING UP DATA
################################################################


# NIMBLE code for the Batch Marking Likelihood
nimLc <- nimbleCode({ 
  ## Priors and constraints
  
  # unmarked params
  lambda ~ dgamma(shape=5,scale=2000) # prior for initial unmarked population, mean = shape * scale, variance = shape * scale^2
  eta    ~ dunif(0, 10) # prior for unmarked recruitment rate
  
  for(t in 1:(To-1)) {
    powEta[t] <- pow(eta, deltaTime[t])
    powPhi[t] <- pow(phi[t], deltaTime[t])
  }
  
  for (t in 1:(To-1)) {
    phi_u[t] <- powPhi[t]
  } #phi
  for (t in 1:To) {
    p_u[t] <- pgt[t,t]
  } #p
  
  # marked params and latent vars
  for (g in 1:G){
    Mgt[g,g] <- Rg[g] # initial marked populations
    
    for (t in 1:(To-1)) {
      phi_m[g,t] <- powPhi[t]
    } #phi
    for (t in 1:(To-1)) {
      p_m[g,t] <- pgt[g,t+1]
    } #t
  } #g
  
  # shared params
  Bphi   ~ dnorm(0,3) # Prior for mean survival
  Bphi_T ~ dnorm(0,3) # Prior for TEMPERATURE covariate coeff
  Bphi_S ~ dnorm(0,3) # Prior for SALINITY covariate coeff
  for(t in 1:(To-1)) {
    phi[t] <- plogis(Bphi + Bphi_S * SALIN[t+1] + Bphi_T * TEMP[t+1])
    #phi[t] <- plogis(Bphi + Bphi_S * SALIN[t+1])
    #phi[t] <- plogis(Bphi)
  }
  
  Bp_T ~ dnorm(0,3) # Prior for TEMPERATURE covariate coeff
  Bp_S ~ dnorm(0,3) # Prior for SALINITY covariate coeff
  Bp_L ~ dnorm(0,3) # Prior for individual LENGTH covariate coeff
  for(t in 1:To) {
    Bp_t[t] ~ dnorm(0,3) # Prior for mean detection (marked)
    
    #p[t] <- plogis(Bp + Bp_S * SALIN[t] + Bp_T * TEMP[t])
    for(g in 1:To) { # note that we use 'To' instead of 'G' here since unmarked lengths are measured at g==To
      pgt[g,t] <- plogis(Bp_t[t] + Bp_S * SALIN[t] + Bp_T * TEMP[t] + Bp_L * MEAN_LENGTH[g,t])
    }
  }
  
  ## Likelihood (unmarked)
  # t==1
  Ut[1] ~ dpois(lambda)
  ut[1] ~ dbinom(p_u[1], Ut[1])
  for (t in 2:To) {
    # State process
    Ut[t] <- At[t] + Ct[t]
    At[t] ~ dbinom(phi_u[t-1], Ut[t-1] - Mgt[t-1,t-1]) # unmarked survival
    Ct[t] ~ dpois(Ut[t-1] * powEta[t-1])               # unmarked recruitment
    
    # Observation process
    ut[t] ~ dbinom(p_u[t], Ut[t])
  } # t
  
  ## Likelihood (marked)
  for (g in 1:G) {
    for (t in (g+1):To) {
      # State process
      Mgt[g,t] ~ dbinom(phi_m[g,t-1], Mgt[g,t-1]) # marked survival; Note: Mgg == Rg
      dgt[g,t] <- Mgt[g,t-1] - Mgt[g,t]           # deaths are those that did not survive, implicitly modelled by Xgt
      
      # Observation process
      mgt[g,t] ~ dbinom(p_m[g,t-1], Mgt[g,t])
    } # t
  } # g
})

# Set up the model
# Assign values to the constants in the model, provide data, initialize the random variables, and then create the model.
LcConsts <- list(To = To,
                 G = G,
                 Sg = Sg,
                 lt = lt,
                 max_i = maxI,
                 n_ind = length(i_weight),
                 deltaTime = deltaT)

LcData <- list(Rg = Rg,
               mgt = mgt,
               TEMP  = temperature,
               MEAN_LENGTH = meanLengths,
               ut = ut,
               SALIN = salinity)

# Initial values
Ut_init = rep(NA, times = To)
Ut_init[1] = 5000
Mgt_init = mgt * NA
Mgt_init[,1] = 5 * mgt[,1]
LcInits <- list(Mgt = Mgt_init,
                Ut = Ut_init, 
                At = Ut_init,
                Ct = 0 * Ut_init,
                phi = rep(1, times=(To-1)), 
                pgt = matrix(ut[1]/Ut_init[1], nrow=To, ncol=To),
                Bp_t = rep(boot::logit(ut[1]/Ut_init[1]), times=To),
                Bp_T = -0.05,
                Bp_S = 0,  
                Bp_L = 0,
                Bphi = 5,
                Bphi_T = -0.05,
                Bphi_S = 0,
                lambda = 5000,
                eta = 0.5,
                p_u = rep(0.05, times=To), 
                phi_u = rep(0.85, times=(To-1)), 
                phi_m = matrix(0.85, nrow=G, ncol=(To-1)),
                TEMP = temperature.inits,
                SALIN = salinity.inits)

# Parameters monitored
parameters <- c("phi", "p_u", "p_m", "lambda", "eta", "Bp_t", "Bp_T", "Bp_S", "Bp_L", "Bphi", "Bphi_T", "Bphi_S")
parameters2 <- c("phi", "pgt", "Bp_t", "Bp_T", "Bp_S", "Bp_L", "Bphi", "Bphi_T", "Bphi_S", "lambda", "eta", "Ut", "At", "Ct", "Mgt", "dgt") # need to monitor all latent states for WAIC calculation

# Create the model
Lc_mod <- nimbleModel(code = nimLc, name = "Lc", data = LcData, constants = LcConsts, inits = LcInits)

# Compile the model in C++
CLc_mod <- compileNimble(Lc_mod)

# Run the MCMC algorithm provided by NIMBLE
LcMCMC  <- buildMCMC(CLc_mod, monitors = parameters, monitors2 = parameters2, enableWAIC = TRUE)
CLcMCMC <- compileNimble(LcMCMC)

tictoc::tic()
LcMCMC_samples <- runMCMC(CLcMCMC, nburnin = 3000000, niter = 3400000, thin = 4, WAIC=TRUE, samplesAsCodaMCMC = T) # 604.99 sec elapsed
tictoc::toc() 

LcMCMC_samples$WAIC

saveRDS(LcMCMC_samples, file = "02 - Case Study 2/MCMC_sandlance_final.RDS")

library(ggmcmc)
ggmcmc::ggs_compare_partial(ggmcmc::ggs(LcMCMC_samples$samples)) + ggplot2::facet_wrap(facets = "Parameter", ncol = 10, scales = "free")
ggmcmc::ggs_traceplot(ggmcmc::ggs(LcMCMC_samples$samples)) + ggplot2::facet_wrap(facets = "Parameter", ncol = 10, scales = "free")

# Parameter Estimates and Credible Intervals, used in producing Table 6 and Table 7:
# colMeans(LcMCMC_samples$samples)
# > colMeans(LcMCMC_samples$samples)
# Bp_L          Bp_S          Bp_T          Bp_t[1]       Bp_t[2]       Bp_t[3]       Bp_t[4] 
# 7.852067e-01 -2.675314e-01  1.935835e+00 -2.008580e+00 -2.695935e+00 -2.152638e+00 -3.117269e+00 
# Bp_t[5]       Bp_t[6]       Bp_t[7]       Bphi          Bphi_S        Bphi_T        eta 
# -8.257290e-01 -2.677747e+00 -2.596188e-01 1.512383e+00  8.765545e-01 -1.827574e+00  2.410084e-01 
# lambda        p_m[1, 1]     p_m[2, 1]     p_m[3, 1]     p_m[4, 1]     p_m[5, 1]     p_m[6, 1] 
# 3.863345e+04  6.948885e-02  4.467023e-02  8.699730e-02  8.699730e-02  8.699730e-02  8.699730e-02 
# p_m[1, 2]     p_m[2, 2]     p_m[3, 2]     p_m[4, 2]     p_m[5, 2]     p_m[6, 2]     p_m[1, 3] 
# 1.048168e-01  1.048168e-01  6.875731e-02  1.048168e-01  1.048168e-01  1.048168e-01  3.505890e-02 
# p_m[2, 3]     p_m[3, 3]     p_m[4, 3]     p_m[5, 3]     p_m[6, 3]     p_m[1, 4]     p_m[2, 4] 
# 2.300245e-02  4.267141e-02  5.990615e-02  4.267141e-02  4.267141e-02  2.540383e-01  1.839709e-01 
# p_m[3, 4]     p_m[4, 4]     p_m[5, 4]     p_m[6, 4]     p_m[1, 5]     p_m[2, 5]     p_m[3, 5] 
# 1.902248e-01  1.957731e-01  2.235774e-01  1.893318e-01  8.275474e-02  8.386302e-02  5.207315e-02 
# p_m[4, 5]     p_m[5, 5]     p_m[6, 5]     p_m[1, 6]     p_m[2, 6]     p_m[3, 6]     p_m[4, 6] 
# 6.833244e-02  6.753938e-02  1.353180e-01  9.249502e-01  9.132327e-01  9.205579e-01  9.137305e-01 
# p_m[5, 6]     p_m[6, 6]     p_u[1]        p_u[2]        p_u[3]        p_u[4]        p_u[5] 
# 9.227016e-01  9.137305e-01  8.657924e-03  4.467023e-02  6.875731e-02  5.990615e-02  2.235774e-01 
# p_u[6]        p_u[7]        phi[1]        phi[2]        phi[3]        phi[4]        phi[5] 
# 1.353180e-01  9.205846e-01  6.341988e-01  8.175234e-01  8.175234e-01  8.897810e-01  8.353523e-01 
# phi[6] 
# 2.108125e-01

# relevant p_m's: p_m[g,t] with g <= t
# p_m[1, 1]    p_m[1, 2]    p_m[1, 3]    p_m[1, 4]    p_m[1, 5]    p_m[1, 6]
# 6.948885e-02 1.048168e-01 3.505890e-02 2.540383e-01 8.275474e-02 9.249502e-01
#              p_m[2, 2]    p_m[2, 3]    p_m[2, 4]    p_m[2, 5]    p_m[2, 6]    
#              1.048168e-01 2.300245e-02 1.839709e-01 8.386302e-02 9.132327e-01
#                           p_m[3, 3]    p_m[3, 4]    p_m[3, 5]    p_m[3, 6] 
#                           4.267141e-02 1.902248e-01 5.207315e-02 9.205579e-01
#                                        p_m[4, 4]    p_m[4, 5]    p_m[4, 6]
#                                        1.957731e-01 6.833244e-02 9.137305e-01
#                                                     p_m[5, 5]    p_m[5, 6] 
#                                                     6.753938e-02 9.227016e-01
#                                                                  p_m[6, 6]
#                                                                  9.137305e-01

# apply(LcMCMC_samples$samples, 2, function(x) quantile(x, 0.025))
# > apply(LcMCMC_samples$samples, 2, function(x) quantile(x, 0.025))
# Bp_L          Bp_S          Bp_T          Bp_t[1]       Bp_t[2]       Bp_t[3]       Bp_t[4] 
# 5.533447e-01 -9.727711e-01  6.372626e-01  -3.041178e+00 -3.457987e+00 -2.447771e+00 -3.372436e+00 
# Bp_t[5]       Bp_t[6]       Bp_t[7]       Bphi          Bphi_S        Bphi_T        eta 
# -1.269121e+00 -2.982714e+00 -1.796454e+00 1.188940e+00  2.925757e-01 -2.221011e+00  1.889889e-01 
# lambda        p_m[1, 1]     p_m[2, 1]     p_m[3, 1]     p_m[4, 1]     p_m[5, 1]     p_m[6, 1] 
# 2.952096e+04  5.101143e-02  3.542638e-02  6.094375e-02  6.094375e-02  6.094375e-02  6.094375e-02 
# p_m[1, 2]     p_m[2, 2]     p_m[3, 2]     p_m[4, 2]     p_m[5, 2]     p_m[6, 2]     p_m[1, 3] 
# 7.960167e-02  7.960167e-02  5.615759e-02  7.960167e-02  7.960167e-02  7.960167e-02  2.732250e-02 
# p_m[2, 3]     p_m[3, 3]     p_m[4, 3]     p_m[5, 3]     p_m[6, 3]     p_m[1, 4]     p_m[2, 4] 
# 1.697071e-02  3.316811e-02  4.518292e-02  3.316811e-02  3.316811e-02  1.899829e-01  1.400329e-01 
# p_m[3, 4]     p_m[4, 4]     p_m[5, 4]     p_m[6, 4]     p_m[1, 5]     p_m[2, 5]     p_m[3, 5] 
# 1.446262e-01  1.486505e-01  1.686590e-01  1.439728e-01  6.312507e-02  6.395077e-02  3.771913e-02 
# p_m[4, 5]     p_m[5, 5]     p_m[6, 5]     p_m[1, 6]     p_m[2, 6]     p_m[3, 6]     p_m[4, 6] 
# 5.173191e-02  5.104495e-02  1.006403e-01  5.338895e-01  4.888910e-01  5.168183e-01  4.906146e-01 
# p_m[5, 6]     p_m[6, 6]     p_u[1]        p_u[2]        p_u[3]        p_u[4]        p_u[5] 
# 5.251364e-01  4.906146e-01  6.621778e-03  3.542638e-02  5.615759e-02  4.518292e-02  1.686590e-01 
# p_u[6]        p_u[7]        phi[1]        phi[2]        phi[3]        phi[4]        phi[5] 
# 1.006403e-01  5.169501e-01  5.242731e-01  7.665514e-01  7.665514e-01  8.465820e-01  7.799146e-01 
# phi[6] 
# 1.783220e-01


# relevant p_m's: p_m[g,t] with g <= t
# p_m[1, 1]    p_m[1, 2]    p_m[1, 3]    p_m[1, 4]    p_m[1, 5]    p_m[1, 6]
# 5.101143e-02 7.960167e-02 2.732250e-02 1.899829e-01 6.312507e-02 5.338895e-01
#              p_m[2, 2]    p_m[2, 3]    p_m[2, 4]    p_m[2, 5]    p_m[2, 6]    
#              7.960167e-02 1.697071e-02 1.400329e-01 6.395077e-02 4.888910e-01
#                           p_m[3, 3]    p_m[3, 4]    p_m[3, 5]    p_m[3, 6] 
#                           3.316811e-02 1.446262e-01 3.771913e-02 5.168183e-01
#                                        p_m[4, 4]    p_m[4, 5]    p_m[4, 6]
#                                        1.486505e-01 5.173191e-02 4.906146e-01
#                                                     p_m[5, 5]    p_m[5, 6] 
#                                                     5.104495e-02 5.251364e-01
#                                                                  p_m[6, 6]
#                                                                  4.906146e-01

# apply(LcMCMC_samples$samples, 2, function(x) quantile(x, 0.975))
# > apply(LcMCMC_samples$samples, 2, function(x) quantile(x, 0.975))
# Bp_L          Bp_S          Bp_T          Bp_t[1]       Bp_t[2]       Bp_t[3]       Bp_t[4] 
# 1.011052e+00  4.048435e-01  3.006847e+00  -1.073506e+00 -1.964950e+00 -1.895726e+00 -2.895379e+00 
# Bp_t[5]       Bp_t[6]       Bp_t[7]       Bphi          Bphi_S        Bphi_T        eta 
# -4.111506e-01 -2.371474e+00 1.036988e+00  2.067396e+00  1.397807e+00 -1.423120e+00  3.139783e-01 
# lambda        p_m[1, 1]     p_m[2, 1]     p_m[3, 1]     p_m[4, 1]     p_m[5, 1]     p_m[6, 1] 
# 4.838159e+04  8.889539e-02  5.457012e-02  1.157057e-01  1.157057e-01  1.157057e-01  1.157057e-01 
# p_m[1, 2]     p_m[2, 2]     p_m[3, 2]     p_m[4, 2]     p_m[5, 2]     p_m[6, 2]     p_m[1, 3] 
# 1.305930e-01  1.305930e-01  8.123056e-02  1.305930e-01  1.305930e-01  1.305930e-01  4.380566e-02 
# p_m[2, 3]     p_m[3, 3]     p_m[4, 3]     p_m[5, 3]     p_m[6, 3]     p_m[1, 4]     p_m[2, 4] 
# 3.070612e-02  5.238248e-02  7.377217e-02  5.238248e-02  5.238248e-02  3.123780e-01  2.248866e-01 
# p_m[3, 4]     p_m[4, 4]     p_m[5, 4]     p_m[6, 4]     p_m[1, 5]     p_m[2, 5]     p_m[3, 5] 
# 2.324718e-01  2.391788e-01  2.734736e-01  2.313880e-01  1.030349e-01  1.042911e-01  6.927560e-02 
# p_m[4, 5]     p_m[5, 5]     p_m[6, 5]     p_m[1, 6]     p_m[2, 6]     p_m[3, 6]     p_m[4, 6] 
# 8.684182e-02  8.599103e-02  1.672115e-01  9.938335e-01  9.924823e-01  9.933345e-01  9.925432e-01 
# p_m[5, 6]     p_m[6, 6]     p_u[1]        p_u[2]        p_u[3]        p_u[4]        p_u[5] 
# 9.935820e-01  9.925432e-01  1.141723e-02  5.457012e-02  8.123056e-02  7.377217e-02  2.734736e-01 
# p_u[6]        p_u[7]        phi[1]        phi[2]        phi[3]        phi[4]        phi[5] 
# 1.672115e-01  9.933379e-01  7.557103e-01  8.876936e-01  8.876936e-01  9.413585e-01  9.053318e-01 
# phi[6] 
# 2.876456e-01


# relevant p_m's: p_m[g,t] with g <= t
# p_m[1, 1]    p_m[1, 2]    p_m[1, 3]    p_m[1, 4]    p_m[1, 5]    p_m[1, 6]
# 8.889539e-02 1.305930e-01 4.380566e-02 3.123780e-01 1.030349e-01 9.938335e-01
#              p_m[2, 2]    p_m[2, 3]    p_m[2, 4]    p_m[2, 5]    p_m[2, 6]    
#              1.305930e-01 3.070612e-02 2.248866e-01 1.042911e-01 9.924823e-01
#                           p_m[3, 3]    p_m[3, 4]    p_m[3, 5]    p_m[3, 6] 
#                           5.238248e-02 2.324718e-01 6.927560e-02 9.933345e-01
#                                        p_m[4, 4]    p_m[4, 5]    p_m[4, 6]
#                                        2.391788e-01 8.684182e-02 9.925432e-01
#                                                     p_m[5, 5]    p_m[5, 6] 
#                                                     8.599103e-02 9.935820e-01
#                                                                  p_m[6, 6]
#                                                                  9.925432e-01

# population size posterior estimates
df = NULL
Nhat_df = NULL
for(j in 1:7) {
  nam <- paste("Ut[", j, "]", sep = "")
  assign(nam, LcMCMC_samples$samples2[, nam])
  df[[nam]] = c(get(nam))
  Nhat_df = rbind(Nhat_df, data.frame(time=j, Nmin=min(df[[nam]]), 
                                      Nmax=max(df[[nam]]), 
                                      Nmedian=median(df[[nam]]), 
                                      Nupper=quantile(df[[nam]], probs = 0.975),
                                      Nlower=quantile(df[[nam]], probs = 0.025),
                                      ut=ut[j]
  ))
}

timeLabels = c(0, cumsum(deltaT[-6])) + 1 # first point is labelled 'Week 1'
Nhat_df$timeWeeks = c(timeLabels, NA)

library(ggplot2)
ggplot(data=Nhat_df) + 
  geom_ribbon(aes(x=timeWeeks, ymin=Nlower, ymax=Nupper, fill="Estimated Unmarked"), colour="blue", fill="blue", alpha=0.1, linewidth=1.15) + # 95% Credible Intervals
  geom_line(aes(x=timeWeeks, y=Nmedian, color="Estimated Unmarked"), linewidth=1.25) +
  geom_point(aes(y=ut, x=timeWeeks, color="Observed Unmarked"), size=1.5) +
  geom_line(aes(y=ut, x=timeWeeks), colour='red', linewidth=1.25, alpha=0.2) +
  scale_x_continuous(labels=1:8, breaks=1:8) +
  xlab("Week") +
  ylab("Population Size") +
  scale_color_manual(
    "Legend", 
    guide = "legend",
    values=c("Estimated Unmarked"="blue", 
             "Observed Unmarked"="red")
  ) +
  scale_fill_manual(
    "Legend", 
    guide = "legend",
    values=c("Estimated Unmarked"="blue", 
             "Observed Unmarked"="red")
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 1)) -> plot.N

ggplot(data=Nhat_df) + 
  geom_errorbar(aes(x=timeWeeks, ymin=Nlower, ymax=Nupper), colour="blue", alpha=0.8, linewidth=1.15) + # 95% Credible Intervals
  geom_point(aes(x=timeWeeks, y=Nmedian, color="Estimated Unmarked"), size=1.5) +
  geom_point(aes(y=ut, x=timeWeeks, color="Observed Unmarked"), size=1.5) +
  scale_x_continuous(labels=1:8, breaks=1:8) +
  xlab("Week") +
  ylab("Population Size") +
  scale_color_manual(
    "Legend", 
    guide = "legend",
    values=c("Estimated Unmarked"="blue", 
             "Observed Unmarked"="red")
  ) +
  scale_fill_manual(
    "Legend", 
    guide = "legend",
    values=c("Estimated Unmarked"="blue", 
             "Observed Unmarked"="red")
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 1)) -> plot.N2


# Add the observed population to the estimated unobserved to get the total estimated population size
library(dplyr)
obs_df = data.frame(N_obs=colSums(mgt, na.rm = TRUE) + ut, time = 1:7)
Ntot_df <- Nhat_df %>%
  inner_join(obs_df, by = "time") %>%
  mutate(
    Nmin = Nmin + N_obs,
    Nmax = Nmax + N_obs,
    Nmedian = Nmedian + N_obs,
    Nupper = Nupper + N_obs,
    Nlower = Nlower + N_obs
  ) %>%
  select(-c(N_obs))

# Produce Figure 3
ggplot(data=Ntot_df) + 
  geom_errorbar(aes(x=timeWeeks, ymin=Nlower, ymax=Nupper), colour="blue", alpha=0.4, linewidth=1.15) + # 95% Credible Intervals
  geom_point(aes(x=timeWeeks, y=Nmedian, color="Total Population"), size=2) +
  geom_point(aes(y=ut, x=timeWeeks, color="Observed Unmarked"), size=2) +
  scale_x_continuous(labels=1:8, breaks=1:8) +
  xlab("Week") +
  ylab("Population Size") +
  scale_color_manual(
    "Legend", 
    guide = "legend",
    values=c("Total Population"="blue", 
             "Observed Unmarked"="red")
  ) +
  scale_fill_manual(
    "Legend", 
    guide = "legend",
    values=c("Total Population"="blue", 
             "Observed Unmarked"="red")
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.8, 1)) -> plot.N3
