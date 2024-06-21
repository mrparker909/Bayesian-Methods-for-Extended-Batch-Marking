library(nimble)

######################################################
#' DEFINITIONS:
#' 
#' OBSERVED DATA:
#' To  the number of capture–recapture occasions.
#' G   the number of batch-marked release groups; G ≤ T.
#' Sg  the number of individuals sampled at sampling occasion g; g = 1, 2, . . . , G.
#' Rg  the number of individuals marked and released at sampling occasion g from batch group g; g = 1, 2, . . . , G. We condition on the {Rg} when we form the likelihood for marked individuals.
#' rgt the number of individuals from batch group g recaptured at recapture occasion t; g = 1, 2, . . . , G, t = g + 1, . . . , T.
#' ut  the number of individuals captured at sampling occasion t that were not marked; t = 1, . . . , T.
#' lt the number of individuals lost at sampling occasion t; t = 2, . . . , T. ℓt = ut − Rt.
#'
#' LATENT VARIABLES:
#' Xgt the number of individuals present at occasion t from marked group g; g = 1, 2, . . . , G, t = g + 1, . . . , T.
#' dgt the number of individuals from marked group g that die at sampling occasion t. 
#' Ut  the total unmarked individuals present at occasion t
#' At  the total unmarked individuals which survive from t to t+1
#' Ct  the total unmarked individuals recruited from t to t+1
#' 
#' PARAMETERS:
#' phi    survival
#' p      capture probability
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

rgt = array(0, dim=c(G,To))
for(g in 1:G) {
  # number of individuals from group g recaptured at recapture occasion t
  rgt[g,] = c(0,countsMat[g,-1])
}
# DONE SETTING UP DATA
################################################################


#NIMBLE code for the Batch Marking Likelihood
nimLm <- nimbleCode({ 
  ## Priors and constraints
  
  # marked params and latent vars
  for (g in 1:G){
    Xgt[g,g] <- Rg[g] # initial marked populations
    
    for (t in 1:(To-1)) {
      phi_m[g,t] <- phi[t]
    } #phi
    for (t in 1:(To-1)) {
      p_m[g,t] <- p[t]
    } #p
  } #i
  
  # shared params
  for(t in 1:(To-1)) {
    phi[t] ~ dunif(0, 1) # Prior for mean survival
  }
  for(t in 1:(To-1)) {
    p[t]   ~ dunif(0, 1) # Prior for mean recapture
  }
  
  ## Likelihood (marked)
  for (g in 1:G){
    for (t in (g+1):To) {
      # State process
      Xgt[g,t] ~ dbinom(phi_m[g,t-1], Xgt[g,t-1]) # marked survival; Note: Xgg == Rg
      #dgt[g,t] <- Xgt[g,t-1] - Xgt[g,t]          # deaths are those that did not survive, implicitly modelled by Xgt
      
      # Observation process
      rgt[g,t] ~ dbinom(p_m[g,t-1], Xgt[g,t])
    } # t
  } # g
})

# Set up the model
# Assign values to the constants in the model, provide data, initialize the random variables, and then create the model.
LmConsts <- list(To = To,
                 G = G,
                 Sg = Sg,
                 Rg = Rg,
                 rgt = rgt)

# Initial values
Xgt_init = rgt * NA
Xgt_init[,1] = 5 * rgt[,1]
LmInits <- list(Xgt = Xgt_init, 
                phi = rep(0.7, times=(To-1)), 
                p = rep(0.1, times=(To-1)),
                p_m = matrix(0.1, nrow=G, ncol=(To-1)), 
                phi_m = matrix(0.5, nrow=G, ncol=(To-1)))

# Parameters monitored
parameters <- c("phi", "p")
parameters2 <- c("phi", "p", "Xgt") # need to monitor all latent states for WAIC calculation

# Create the model
Lm_mod <- nimbleModel(code = nimLm, name = "Lm", constants = LmConsts, inits = LmInits)

# Compile the model in C++
CLm_mod <- compileNimble(Lm_mod)

# Run the MCMC algorithm provided by NIMBLE
tictoc::tic()
box <- list( list(c("phi", "p"), c(0, 1))) ## Constraints for the parameters
LmMCMC <- buildMCMC(CLm_mod, boxConstraints = box, monitors = parameters, monitors2 = parameters2, enableWAIC = TRUE)
CLmMCMC <- compileNimble(LmMCMC)

LmMCMC_samples <- runMCMC(CLmMCMC, nburnin = 50000, niter = 100000, WAIC=TRUE, samplesAsCodaMCMC = T)
LmMCMC_samples$WAIC
tictoc::toc()

saveRDS(LmMCMC_samples, file = "CaseStudy-Cowen2017/Marked Likelihood/Lm-results-phi_t_p_t.RDS")

# ggmcmc::ggs_compare_partial(ggs(LmMCMC_samples$samples)) + facet_wrap(ncol = 4, facets = "Parameter", scales="free")