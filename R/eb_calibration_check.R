#' Goodness-of-Fit Calibration Check for Empirical Bayes Models
#'
#' Performs a simulation-based goodness-of-fit assessment for an empirical
#' Bayes model using the multivariate energy distance as the discrepancy
#' statistic.
#'
#' @param X_obs A numeric matrix or data frame of observed data (n x p).
#'   A numeric vector is coerced to a single-column matrix.
#' @param prior_sampler function(n): draws n latent parameters from G_hat.
#' @param likelihood_sampler function(thetas): draws one observation per
#'   latent parameter, on the natural scale.
#' @param feature_transform function(x): maps natural-scale observations to the
#'   scale on which the energy statistic is computed (variance-stabilizing
#'   transform, or identity). Applied identically to observed and simulated
#'   data. Default: identity.
#' @param B Integer. Number of Monte Carlo replicates. Default: 500.
#' @param seed Optional integer random seed. Default: 1.
#' @param standardize Logical. If TRUE, center/scale each feature column using
#'   statistics computed from the (transformed) observed data, and reuse those
#'   same statistics for the simulated data. Default: TRUE.
#' @param verbose Logical. If TRUE, prints a progress bar. Default: TRUE.
#'
#' @return A list of class "eb_calibration".
#' @export
eb_calibration_check <- function(X_obs,
                                 prior_sampler,
                                 likelihood_sampler,
                                 feature_transform = function(x) x,
                                 B = 500,
                                 seed = 1,
                                 standardize = TRUE,
                                 verbose = TRUE) {

  # ---- Input validation -----------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  if (is.vector(X_obs) && !is.list(X_obs)) {
    X_obs <- matrix(X_obs, ncol = 1L)
  } else {
    X_obs <- as.matrix(X_obs)
  }

  n <- nrow(X_obs)

  if (!is.function(prior_sampler))      stop("`prior_sampler` must be a function.")
  if (!is.function(likelihood_sampler)) stop("`likelihood_sampler` must be a function.")
  if (!is.function(feature_transform))  stop("`feature_transform` must be a function.")
  if (!is.numeric(B) || length(B) != 1L || B < 1L) {
    stop("`B` must be a positive integer.")
  }
  B <- as.integer(B)

  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required. Install it with: install.packages('energy')")
  }

  # ---- Energy distance via the energy package -------------------------------
  energy_distance <- function(X, Y) {
    X  <- as.matrix(X)
    Y  <- as.matrix(Y)
    nx <- nrow(X)
    ny <- nrow(Y)
    XY <- rbind(X, Y)
    e_stat <- energy::eqdist.e(XY, sizes = c(nx, ny))
    e_stat * (nx + ny) / (nx * ny)
  }

  # ---- Feature transform + standardization ----------------------------------
  # Applies feature_transform, then (optionally) centers/scales using
  # reference statistics. When ref_stats is NULL the stats are computed from
  # the data passed in (used once, on the observed data); thereafter the
  # observed-data stats are reused so observed and simulated are on one scale.
  apply_features <- function(x, ref_stats = NULL) {
    Z <- as.matrix(feature_transform(x))
    if (standardize) {
      if (is.null(ref_stats)) {
        mu <- colMeans(Z)
        s  <- apply(Z, 2, stats::sd)
      } else {
        mu <- ref_stats$mu
        s  <- ref_stats$s
      }
      s[s == 0] <- 1
      Z <- sweep(sweep(Z, 2L, mu, "-"), 2L, s, "/")
      attr(Z, "mu") <- mu
      attr(Z, "s")  <- s
    }
    Z
  }

  # Observed features + reference standardization stats (from observed only)
  Z_obs     <- apply_features(X_obs)
  ref_stats <- if (standardize) {
    list(mu = attr(Z_obs, "mu"), s = attr(Z_obs, "s"))
  } else NULL

  # ---- Helper: one synthetic dataset, on the feature scale ------------------
  generate_synthetic <- function() {
    thetas <- prior_sampler(n)
    sim    <- likelihood_sampler(thetas)
    if (is.vector(sim) && !is.list(sim)) sim <- matrix(sim, ncol = 1L)
    apply_features(sim, ref_stats)
  }

  # ---- Main Monte Carlo loop ------------------------------------------------
  T_obs_b <- numeric(B)
  T_ref_b <- numeric(B)

  if (verbose) {
    cat("Running empirical Bayes calibration check (B =", B, "replicates)...\n")
    pb <- utils::txtProgressBar(min = 0L, max = B, style = 3L)
  }

  for (b in seq_len(B)) {
    sim1 <- generate_synthetic()
    sim2 <- generate_synthetic()

    T_ref_b[b] <- energy_distance(sim1, sim2)
    T_obs_b[b] <- energy_distance(Z_obs, sim1)

    if (verbose) utils::setTxtProgressBar(pb, b)
  }

  if (verbose) {
    close(pb)
    cat("\n")
  }

  # ---- Test statistic and p-value -------------------------------------------
  T_obs_median <- stats::median(T_obs_b)
  p_value      <- (1 + sum(T_ref_b >= T_obs_median)) / (B + 1)

  structure(
    list(
      p_value   = p_value,
      T_obs     = T_obs_median,
      T_obs_all = T_obs_b,
      T_ref     = T_ref_b,
      B         = B
    ),
    class = "eb_calibration"
  )
}


# ---- S3 methods -------------------------------------------------------------

#' @export
print.eb_calibration <- function(x, ...) {
  cat("Empirical Bayes Goodness-of-Fit Calibration Check\n")
  cat("--------------------------------------------------\n")
  cat(sprintf("  Replicates (B)       : %d\n", x$B))
  cat(sprintf("  Median T_obs         : %.4f\n", x$T_obs))
  cat(sprintf("  Monte Carlo p-value  : %.4f\n", x$p_value))
  cat("\n")
  if (x$p_value < 0.05) {
    cat("  >> Evidence of poor fit (p < 0.05).\n")
  } else {
    cat("  >> No strong evidence against the fitted model (p >= 0.05).\n")
  }
  invisible(x)
}


#' @export
plot.eb_calibration <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(1L, 2L), mar = c(4, 4, 3, 1))

  # --- Panel 1: Reference distribution vs T_obs --------------------------------
  hist(x$T_ref,
       breaks  = 30L,
       col     = "steelblue",
       border  = "white",
       xlab    = "Energy Distance",
       main    = "Reference Distribution vs T_obs",
       probability = TRUE,
       ...)
  abline(v   = x$T_obs,
         col = "firebrick",
         lwd = 2L,
         lty = 2L)
  legend("topright",
         legend = c(
           expression(T^{ref}),
           bquote(tilde(T)^{obs} == .(round(x$T_obs, 3)))
         ),
         col    = c("steelblue", "firebrick"),
         lwd    = c(8L, 2L),
         lty    = c(1L, 2L),
         bty    = "n")

  # --- Panel 2: Per-replicate observed vs reference discrepancies --------------
  plot(x$T_ref, x$T_obs_all,
       pch  = 16L,
       cex  = 0.5,
       col  = adjustcolor("steelblue", alpha.f = 0.5),
       xlab = expression(T[b]^{ref}),
       ylab = expression(T[b]^{obs}),
       main = "Per-Replicate Discrepancies")
  abline(0, 1, col = "grey40", lty = 2L)
  abline(h = x$T_obs, col = "firebrick", lwd = 1.5, lty = 3L)
  mtext(sprintf("p = %.4f", x$p_value),
        side = 3L, line = -1.5, adj = 0.98, cex = 0.85, col = "firebrick")

  invisible(x)
}
