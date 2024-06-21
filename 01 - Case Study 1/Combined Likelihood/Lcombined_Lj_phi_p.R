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
# data from Cowen 2017
counts = read.csv(file = "CaseStudy-Cowen2017/Cowen2017_data.csv")
countsMat = as.matrix(counts)
To = 11
G  = 10
Rg = c(countsMat[,1],0)
ut = c(306,187,139,147,89,76,20,52,63,55,38) 
lt = c(26,48,24,21,9,11,6,2,9,0,38)
Sg = c(306,219,189,207,111,96,30,68,83,81,50)

mgt = array(0, dim=c(G,To))
for(g in 1:G) {
  # number of individuals from group g recaptured at recapture occasion t
  mgt[g,] = c(0,countsMat[g,-1])
}
# DONE SETTING UP DATA
################################################################


#NIMBLE code for the Batch Marking Likelihood
nimLc <- nimbleCode({ 
  ## Priors and constraints
  
  # unmarked params
  lambda ~ dgamma(shape=10,scale=100) # initial unmarked population
  eta    ~ dunif(0, 10)               # prior for unmarked recruitment rate
  for (t in 1:(To-1)) {
    phi_u[t] <- phi
  } #phi
  for (t in 1:To) {
    p_u[t] <- p
  } #p
  
  # marked params and latent vars
  for (g in 1:G){
    Mgt[g,g] <- Rg[g] # initial marked populations
    
    for (t in 1:(To-1)) {
      phi_m[g,t] <- phi
    } #phi
    for (t in 1:(To-1)) {
      p_m[g,t] <- p
    } #p
  } #i
  
  # shared params
  phi ~ dunif(0, 1) # Prior for mean survival
  p   ~ dunif(0, 1) # Prior for mean recapture
  
  ## Likelihood (unmarked)
  # t==1
  Ut[1] ~ dpois(lambda)
  ut[1] ~ dbinom(p_u[1], Ut[1])
  for (t in 2:To) {
    # State process
    Ut[t] <- At[t] + Ct[t]
    At[t] ~ dbinom(phi_u[t-1], Ut[t-1] - Mgt[t-1,t-1]) # unmarked survival
    Ct[t] ~ dpois(Ut[t-1] * eta)                       # unmarked recruitment
    
    # Observation process
    ut[t] ~ dbinom(p_u[t], Ut[t])
  } # t
  
  ## Likelihood (marked)
  for (g in 1:G){
    for (t in (g+1):To) {
      # State process
      Mgt[g,t] ~ dbinom(phi_m[g,t-1], Mgt[g,t-1]) # marked survival; Note: Mgg == Rg
      dgt[g,t] <- Mgt[g,t-1] - Mgt[g,t]           # deaths are those that did not survive, implicitly modelled by Mgt
      
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
                 Rg = Rg,
                 mgt = mgt,
                 lt = lt,
                 ut = ut)

# Initial values
Ut_init = rep(NA, times = To)
Ut_init[1] = 5 * Rg[1]
Mgt_init = mgt * NA
Mgt_init[,1] = 5 * mgt[,1]
LcInits <- list(Mgt = Mgt_init, 
                #dgt = 0 * mgt,
                Ut = Ut_init, 
                At = Ut_init,
                Ct = 0 * Ut_init,
                phi = 0.7, 
                p = 0.1, 
                lambda = 5 * Rg[1],
                eta = 1,
                p_u = rep(0.1, times=To), 
                phi_u = rep(0.5, times=(To-1)),
                p_m = matrix(0.1, nrow=G, ncol=(To-1)), 
                phi_m = matrix(0.5, nrow=G, ncol=(To-1)))

# Parameters monitored
parameters <- c("phi", "p", "lambda", "eta")
parameters2 <- c("phi", "p", "lambda", "eta", "Ut", "At", "Ct", "Mgt", "dgt") # need to monitor all latent states for WAIC calculation

# Create the model
Lc_mod <- nimbleModel(code = nimLc, name = "Lc", constants = LcConsts, inits = LcInits)

# Compile the model in C++
CLc_mod <- compileNimble(Lc_mod)

# Run the MCMC algorithm provided by NIMBLE
tictoc::tic()
box <- list( list(c("phi", "p"), c(0, 1))) ## Constraints for the parameters
LcMCMC <- buildMCMC(CLc_mod, boxConstraints = box, monitors = parameters, monitors2 = parameters2, enableWAIC = TRUE)
CLcMCMC <- compileNimble(LcMCMC)

LcMCMC_samples <- runMCMC(CLcMCMC, nburnin = 100000, niter = 200000, WAIC=TRUE, samplesAsCodaMCMC = T)
LcMCMC_samples$WAIC
tictoc::toc()

saveRDS(LcMCMC_samples, file = "CaseStudy-Cowen2017/Combined Likelihood/Lj-results-phi_p_100k.RDS")

# library(ggmcmc)
# saved as 850x300 png
# ggmcmc::ggs_compare_partial(ggs(LcMCMC_samples$samples)) + facet_wrap(ncol = 4, facets = "Parameter", scales="free")