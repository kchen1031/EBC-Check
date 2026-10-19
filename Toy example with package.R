# Install and load
remotes::install_github("kchen1031/EBC-Check")
library(ebcalibration)


# Poisson-Gamma data
set.seed(123)
n <- 500
theta_true <- rgamma(n, shape = 3, rate = 1)
Y_obs <- rpois(n, lambda = theta_true)

likelihood_sampler <- function(thetas) rpois(length(thetas), lambda = pmax(thetas, 0))

# The test refits the prior on every reference dataset, so it needs the
# g-modeling routine itself, not just a sampler from one fitted G_hat.
# fit_prior(x) must return a function(n) drawing from the prior fitted to x.
# Everything fit_prior estimates must come from its argument x, never from
# Y_obs. Unit-specific values (standard errors, exposures, ...) would be
# captured from the enclosing scope and used by position; there are none here.
#
# B = 199 makes alpha * (B + 1) an integer at the 1%, 5% and 10% levels.


# Case 1: Correctly specified prior family (Gamma, fit by moments)
fit_gamma <- function(x) {
  m <- mean(x)
  v <- max(var(x) - m, 1e-6)            # Var(theta) = Var(X) - E(X)
  shape <- m^2 / v; rate <- m / v
  function(n) rgamma(n, shape = shape, rate = rate)
}

result_good <- eb_calibration_check(
  X_obs              = Y_obs,
  fit_prior          = fit_gamma,
  likelihood_sampler = likelihood_sampler,
  feature_transform  = sqrt,
  B = 199, K = 5, seed = 1
)
print(result_good)
plot(result_good)


# Case 2: Misspecified prior family (a point mass, fit by the sample mean)
fit_point <- function(x) {
  m <- mean(x)
  function(n) rep(m, n)
}

result_bad <- eb_calibration_check(
  X_obs              = Y_obs,
  fit_prior          = fit_point,
  likelihood_sampler = likelihood_sampler,
  feature_transform  = sqrt,
  B = 199, K = 5, seed = 1
)
print(result_bad)
plot(result_bad)

# If any refit errors, result$p_value counts it as exceeding T_obs
# (conservative); result$p_value_dropped is the version that ignores it.
c(conservative = result_bad$p_value, dropped = result_bad$p_value_dropped)
