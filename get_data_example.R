get_data <- function(paper, dataset) {
  
  key <- paste0(tolower(trimws(paper)), "_", tolower(trimws(dataset)))
  
  # ------------------------------------------------------------------
  # 4. Narasimhan et al. 2020 - Surgery (deconvolveR package)
  #    Model: X_i | theta_i ~ Binomial(n_i, theta_i), theta_i ~ G
  #    Transform: asin(sqrt(x / n)). Arcsine-sqrt is the binomial VST; it also
  #    puts counts from different trial sizes n_i on a comparable scale.
  # ------------------------------------------------------------------
  if (key == "narasimhan_surgery") {
    
    if (!requireNamespace("deconvolveR", quietly = TRUE))
      stop("Please install the 'deconvolveR' package.")
    
    data("surg", package = "deconvolveR", envir = environment())
    n_obs <- surg$n   # nodes removed (trials)
    x_obs <- surg$s   # malignant nodes (successes)
    
    tau <- seq(0.01, 0.99, length.out = 200)
    
    # Binomial family needs a two-column matrix: cbind(trials, successes)
    X_mat <- cbind(surg$n, surg$s)
    fit   <- deconvolveR::deconv(tau = tau, X = X_mat, family = "Binomial")
    
    G_hat <- list(support = tau, weights = fit$stats[, "g"])
    
    return(list(
      data  = surg,
      G_hat = G_hat,
      
      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },
      
      likelihood_sampler = function(theta) {
        n <- length(theta)
        stats::rbinom(n, size = rep_len(n_obs, n), prob = theta)
      },
      
      # Arcsine-sqrt binomial VST. Closes over the observed trial sizes n_obs,
      # which match the design points used by likelihood_sampler.
      feature_transform = function(x) {
        n <- rep_len(n_obs, length(x))
        asin(sqrt(x / pmax(n, 1e-12)))
      }
    ))
  }
}