# Install and load
remotes::install_github("kchen1031/EBC-Check")
library(ebcalibration)


# Specify Poisson-Gamma model
set.seed(123)
n <- 500
theta_true <- rgamma(n, shape = 3, rate = 1)
Y_obs <- rpois(n, lambda = theta_true)


# Case 1: Correctly specified prior 
prior_sampler_good      <- function(n) rgamma(n, shape = 3, rate = 1)
likelihood_sampler_good <- function(thetas) rpois(length(thetas), lambda = thetas)

result_good <- eb_calibration_check(
  X_obs              = Y_obs,
  prior_sampler      = prior_sampler_good,
  likelihood_sampler = likelihood_sampler_good,
  B                  = 500,
  seed               = 1
)

print(result_good)
plot(result_good)


# Case 2: Misspecified prior
prior_sampler_bad <- function(n) rnorm(n, mean = 3, sd = 0.5)
likelihood_sampler_bad <- function(thetas) {
  rpois(length(thetas), lambda = pmax(thetas, 0.01))
}

result_bad <- eb_calibration_check(
  X_obs              = Y_obs,
  prior_sampler      = prior_sampler_bad,
  likelihood_sampler = likelihood_sampler_bad,
  B                  = 500,
  seed               = 1
)

print(result_bad)
plot(result_bad)
