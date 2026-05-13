#' Goodness-of-Fit Calibration Check for Empirical Bayes Models
#'
#' Performs a simulation-based goodness-of-fit assessment for an empirical
#' Bayes model using the multivariate energy distance as the discrepancy
#' statistic. See: Efron (2010), "Large-Scale Inference", and the procedure
#' described in the accompanying methodology document.
#'
#' @param X_obs A numeric matrix or data frame of observed data (n x p),
#'   where rows are observations and columns are variables. A numeric vector
#'   is coerced to a single-column matrix.
#' @param prior_sampler A function with signature \code{function(n)} that
#'   draws \code{n} latent parameters from the estimated prior \eqn{\hat{G}}.
#'   Returns a vector (univariate) or matrix with \code{n} rows (multivariate).
#' @param likelihood_sampler A function with signature
#'   \code{function(thetas)} that draws one observation per latent parameter.
#'   \code{thetas} is the output of \code{prior_sampler}. Returns a vector or
#'   matrix with the same number of rows as \code{thetas}.
#' @param B Integer. Number of Monte Carlo replicates. Default: \code{999}.
#' @param seed Optional integer random seed for reproducibility. Default: \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints a progress bar. Default: \code{TRUE}.
#'
#' @return A list of class \code{"eb_calibration"} containing:
#'   \describe{
#'     \item{\code{p_value}}{Monte Carlo p-value. Small values indicate
#'       poor fit.}
#'     \item{\code{T_obs}}{Median observed discrepancy \eqn{\tilde{T}^{obs}}.}
#'     \item{\code{T_obs_all}}{Vector of length \code{B} of per-replicate
#'       observed discrepancies \eqn{T^{obs}_b}.}
#'     \item{\code{T_ref}}{Vector of length \code{B} of reference
#'       discrepancies \eqn{T^{ref}_b}.}
#'     \item{\code{B}}{Number of replicates used.}
#'   }
#'
#' @details
#' The procedure works as follows:
#'
#' For each replicate \eqn{b = 1, \ldots, B}:
#' \enumerate{
#'   \item Draw latent parameters \eqn{\theta_i^{(b)} \sim \hat{G}}, for
#'     \eqn{i = 1, \ldots, n}.
#'   \item Generate two independent synthetic datasets
#'     \eqn{X^{(b)}_{sim,1}} and \eqn{X^{(b)}_{sim,2}} from the fitted model.
#'   \item Compute \eqn{T^{ref}_b = d(X^{(b)}_{sim,1}, X^{(b)}_{sim,2})}.
#'   \item Compute \eqn{T^{obs}_b = d(X^{obs}, X^{(b)}_{sim,1})}.
#' }
#'
#' The test statistic is \eqn{\tilde{T}^{obs} = \text{median}_b(T^{obs}_b)}
#' and the Monte Carlo p-value is:
#' \deqn{p = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\{T^{ref}_b \geq \tilde{T}^{obs}\}}{B + 1}}
#'
#' The energy distance between two samples is computed via
#' \code{eqdist.e()} from the \pkg{energy} package.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#'
#' # Simulate observed data from a Normal-Normal model
#' # True prior: theta_i ~ N(0, 1); likelihood: X_i | theta_i ~ N(theta_i, 1)
#' thetas_true <- rnorm(n, mean = 0, sd = 1)
#' X_obs       <- matrix(rnorm(n, mean = thetas_true, sd = 1), ncol = 1)
#'
#' # Suppose the empirical Bayes estimate of G is N(0, 1) (correctly specified)
#' prior_sampler      <- function(n) rnorm(n, mean = 0, sd = 1)
#' likelihood_sampler <- function(thetas) matrix(rnorm(length(thetas),
#'                                                     mean = thetas, sd = 1),
#'                                               ncol = 1)
#'
#' result <- eb_calibration_check(X_obs, prior_sampler, likelihood_sampler,
#'                                B = 499, seed = 1)
#' print(result)
#' plot(result)
#' @importFrom energy eqdist.e
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
eb_calibration_check <- function(X_obs,
                                 prior_sampler,
                                 likelihood_sampler,
                                 B       = 999L,
                                 seed    = NULL,
                                 verbose = TRUE) {

  # ---- Input validation -------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  if (is.vector(X_obs) && !is.list(X_obs)) {
    X_obs <- matrix(X_obs, ncol = 1L)
  } else {
    X_obs <- as.matrix(X_obs)
  }

  n <- nrow(X_obs)

  if (!is.function(prior_sampler))      stop("`prior_sampler` must be a function.")
  if (!is.function(likelihood_sampler)) stop("`likelihood_sampler` must be a function.")
  if (!is.numeric(B) || length(B) != 1L || B < 1L) {
    stop("`B` must be a positive integer.")
  }
  B <- as.integer(B)

  # ---- Energy distance via the energy package ---------------------------------
  # eqdist.e() computes the E-statistic = d(X,Y) * nx*ny/(nx+ny) for a pooled
  # matrix and a size vector. We recover the plain energy distance by dividing
  # by the factor nx*ny/(nx+ny).
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop("Package 'energy' is required. Install it with: install.packages('energy')")
  }

  energy_distance <- function(X, Y) {
    X  <- as.matrix(X)
    Y  <- as.matrix(Y)
    nx <- nrow(X)
    ny <- nrow(Y)
    XY <- rbind(X, Y)
    # eqdist.e returns the E-statistic = d(X,Y) * nx*ny/(nx+ny)
    e_stat <- eqdist.e(XY, sizes = c(nx, ny))
    e_stat * (nx + ny) / (nx * ny)
  }

  # ---- Helper: generate one synthetic dataset from the fitted model -----------
  generate_synthetic <- function() {
    thetas <- prior_sampler(n)
    sim    <- likelihood_sampler(thetas)
    if (is.vector(sim) && !is.list(sim)) sim <- matrix(sim, ncol = 1L)
    as.matrix(sim)
  }

  # ---- Main Monte Carlo loop --------------------------------------------------
  T_obs_b <- numeric(B)
  T_ref_b <- numeric(B)

  if (verbose) {
    cat("Running empirical Bayes calibration check (B =", B, "replicates)...\n")
    pb <- txtProgressBar(min = 0L, max = B, style = 3L)
  }

  for (b in seq_len(B)) {
    sim1 <- generate_synthetic()
    sim2 <- generate_synthetic()

    T_ref_b[b] <- energy_distance(sim1, sim2)
    T_obs_b[b] <- energy_distance(X_obs, sim1)

    if (verbose) setTxtProgressBar(pb, b)
  }

  if (verbose) {
    close(pb)
    cat("\n")
  }

  # ---- Test statistic and p-value ---------------------------------------------
  T_obs_median <- median(T_obs_b)
  p_value      <- (1 + sum(T_ref_b >= T_obs_median)) / (B + 1)

  # ---- Return -----------------------------------------------------------------
  result <- structure(
    list(
      p_value   = p_value,
      T_obs     = T_obs_median,
      T_obs_all = T_obs_b,
      T_ref     = T_ref_b,
      B         = B
    ),
    class = "eb_calibration"
  )

  result
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
#' @importFrom graphics par hist abline legend mtext
#' @importFrom grDevices adjustcolor
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
