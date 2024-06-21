# install.packages("extBatchMarking")
library("extBatchMarking")

###############################################################################

# Examples for both marked and unmarked

#-------------------------------------------------
# Model 1: 10 phi and 10 probs
#-------------------------------------------------
theta <- c(rep(0, 10), rep(-1, 10))

start.time <- Sys.time()
res <- batchMarkOptim(par=theta, data=WeatherLoach, choiceModel = "model1",
                      method="BFGS", parallel=FALSE, hessian = TRUE, control = list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$phi
res$p
res$ll
res$hessian
res$AIC

#-------------------------------------------------
# Model 2: 10 phis and 1 prob
#-------------------------------------------------
theta <- c(0, rep(-1, 10))

start.time <- Sys.time()
res <- batchMarkOptim(par=theta, data=WeatherLoach, choiceModel = "model2", method="BFGS",
                      parallel=FALSE, hessian = TRUE, control = list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

str(res)
res$phi
res$p
res$ll
res$hessian
res$AIC

#-------------------------------------------------
# Model 3: 1 phi and 10 probs
#-------------------------------------------------
theta <- c(0, rep(-1, 10))

start.time <- Sys.time()
res <- batchMarkOptim(par=theta, data=WeatherLoach, choiceModel = "model3", method="BFGS",
                      parallel=FALSE, hessian = TRUE, control = list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$phi
res$p
res$ll
res$hessian
res$AIC

#-------------------------------------------------
# Model 4: 1 phi and 1 prob
#-------------------------------------------------
theta <- c(0, -1)

start.time <- Sys.time()
res <- batchMarkOptim(par=theta, data=WeatherLoach, choiceModel = "model4", method="BFGS",
                      parallel=FALSE, hessian = TRUE, control = list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$phi
res$p
res$ll
res$hessian
res$AIC

########################################################
# Combined of Marked and Unmarked models

#-------------------------------------------------
# # Model 1: 11 phi and 11 probs
#-------------------------------------------------
thet <- c(rep(0.1, 21),7, -1.5)

start.time <- Sys.time()
res <- batchMarkUnmarkOptim(par=thet, WeatherLoach, Umax=1800, nBins=20, choiceModel="model1", method="BFGS", popSize = "Model-Based",
                            parallel=FALSE, control=list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$phi
res$p
res$lambda
res$gam
res$ll
res$AIC

#-------------------------------------------------
# Model 2: 11 phis and 1 probs
#-------------------------------------------------
thet <- c(rep(0.1, 11), 7, -1.5)

start.time <- Sys.time()
res <- batchMarkUnmarkOptim(par=thet, WeatherLoach, Umax=1800, nBins=20, choiceModel="model2", method="BFGS", popSize = "Model-Based",
                            parallel=FALSE, control=list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$N
res$phi
res$p
res$lambda
res$gam
res$ll
res$AIC

#-------------------------------------------------
# # Model 3: 1 phi and 11 probs
#-------------------------------------------------

thet <- c(0.52, 0.538, -1.50, -1.09, -0.83, -2.09, -2.05, -2.65, -1.71, -1.35, -0.96, -1.91, 7, -1.5)

start.time <- Sys.time()
res <- batchMarkUnmarkOptim(par=thet, WeatherLoach, Umax=1800, nBins=20, choiceModel="model3", method="BFGS", popSize = "Model-Based",
                            parallel=FALSE, control=list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

res$phi
res$p
res$lambda
res$gam
res$ll
res$AIC
res$N

#-------------------------------------------------
# Model 4: 1 phi and 1 prob
#-------------------------------------------------
thet <- c(0.1, 0.1, 7, -1.5)

start.time <- Sys.time()
res <- batchMarkUnmarkOptim(par=thet, WeatherLoach, Umax=1800, nBins=20, choiceModel="model4", method="BFGS", popSize = "Model-Based",
                            parallel=FALSE, control=list(trace = 1))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

str(res)
res$phi
res$p
res$lambda
res$gam
res$ll
res$AIC
res$N