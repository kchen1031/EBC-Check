#' Goodness-of-Fit Calibration Check for Empirical Bayes Models
#'
#' Tests whether the observed data are a plausible realisation of the fitted
#' empirical Bayes data generating process, \eqn{H_0: X_1, \dots, X_n \sim
#' f_{\hat G}}, using the energy distance as the discrepancy and a refitted
#' Monte Carlo reference distribution (proposal Section 4.2, Algorithm 1).
#'
#' One map is applied to every dataset \eqn{Y}: fit the prior to \eqn{Y} with
#' `fit_prior`, draw `K` pseudo-samples of size \eqn{n} from that fit, and
#' return the median of the `K` energy distances between \eqn{Y} and its
#' pseudo-samples. \eqn{T^{obs}} is this map applied to `X_obs`; each
#' \eqn{T^{ref}_b} is the same map applied to a synthetic dataset drawn from
#' the fit to `X_obs`, which refits the prior. The test is asymptotically, not
#' exactly, valid, because the reference datasets come from \eqn{f_{\hat G}}
#' rather than the true law.
#'
#' @section Reading the result:
#' A rejection says the fitted prior family does not reproduce the observed
#' marginal. A flexible prior class absorbs many likelihood misspecifications
#' (wrong noise scale, heavy-tailed errors), and dependence among units leaves
#' the marginal unchanged, so a non-rejection certifies neither the likelihood
#' nor independence (proposal Section 4.2.3).
#'
#' @section Unit-specific nuisance values:
#' Standard errors \eqn{s_i}, exposures, trial counts \eqn{N_i} and similar
#' values \eqn{\nu_i} enter through closures. `fit_prior`,
#' `likelihood_sampler` and `feature_transform` should each capture the
#' vector of \eqn{\nu_i} and use it by position: the \eqn{i}th unit of every
#' dataset, observed or synthetic, has the \eqn{i}th \eqn{\nu_i}. In
#' particular `fit_prior(x)` fits using `x` together with the captured
#' \eqn{\nu_i} (e.g. `ebnm(x, s = s_obs)`). If the prior itself depends on
#' \eqn{\nu_i}, the \eqn{i}th draw of the prior sampler must use the \eqn{i}th
#' \eqn{\nu_i}. Every dataset has `n` units in the original order, so the
#' pairing is preserved throughout.
#'
#' @param X_obs Observed data: a numeric vector, or a matrix / data frame
#'   with one row per unit.
#' @param fit_prior `function(x)`: the g-modeling routine. Given a dataset on
#'   the natural scale (a numeric vector when there is one column, otherwise
#'   a matrix), it fits the model and returns either
#'   \itemize{
#'     \item a prior sampler `function(n)` drawing \eqn{n} values of
#'       \eqn{\theta} from the fitted \eqn{\hat G}, or
#'     \item a list with element `prior_sampler` and, optionally,
#'       `likelihood_sampler`, for models whose likelihood also carries
#'       estimated parameters (e.g. a Weibull shape) that must be refit too.
#'   }
#'   It is called once on `X_obs` and once on each of the `B` synthetic
#'   datasets. It should capture the nuisance values \eqn{\nu_i} (see below)
#'   but must not capture `X_obs` or anything fitted to it: everything it
#'   estimates must come from its argument `x`.
#' @param likelihood_sampler `function(theta)`: draws one observation per
#'   latent parameter, on the natural scale, under the working likelihood.
#'   Used whenever the fit does not supply its own.
#' @param feature_transform `function(x)`: maps natural-scale data to the
#'   scale on which the energy distance is computed (variance-stabilising
#'   transform, augmentation with \eqn{\nu_i}, or identity). Default identity.
#' @param B Number of reference datasets, i.e. number of refits. Default 199,
#'   so that \eqn{\alpha(B + 1)} is an integer at the 1%, 5% and 10% levels.
#' @param K Number of pseudo-samples per dataset over which the median
#'   discrepancy is taken, on both sides. Default 5. `K = 1` gives a single
#'   matched replicate.
#' @param alpha Level at which `print()` reports a rejection. Default 0.05.
#' @param standardize Logical. If `TRUE`, each feature column of a dataset is
#'   centred and scaled by that dataset's own mean and sd, and its
#'   pseudo-samples are scaled by the same constants. Default `TRUE`.
#' @param cores Number of cores for the `B` refits (via
#'   [parallel::mclapply()]; not available on Windows). Results are identical
#'   for any value of `cores` given the same `seed`. Default 1.
#' @param seed Optional integer seed for the test's own random numbers. The
#'   caller's random number stream is restored on exit, so code run after the
#'   call is unaffected. With `seed = NULL` the test draws from, and advances,
#'   the caller's stream. Default 1.
#' @param verbose Logical; print progress. Default `TRUE`.
#'
#' @return An object of class `"eb_calibration"`: a list with
#' \describe{
#'   \item{p_value}{Monte Carlo p-value, counting every failed replicate as
#'     exceeding \eqn{T^{obs}}: \eqn{(1 + \#\{T^{ref}_b \ge T^{obs}\} +
#'     \#\mathrm{failed}) / (B + 1)}. Conservative; equals the usual p-value
#'     when nothing fails.}
#'   \item{p_value_dropped}{The p-value computed on successful replicates
#'     only. Anti-conservative if refits fail more often on extreme datasets.}
#'   \item{T_obs}{Observed statistic (median of `T_obs_k`).}
#'   \item{T_obs_k}{The `K` observed-side discrepancies.}
#'   \item{T_ref}{Reference statistics from successful replicates.}
#'   \item{B, B_eff, K}{Requested and successful replicates, inner draws.}
#'   \item{n_failed}{Replicates in which the refit or simulation errored.}
#'   \item{alpha, n}{Reporting level and sample size.}
#' }
#'
#' @examples
#' \dontrun{
#' # Normal means with known unit-specific s_i, N(mu, tau^2) prior by moments.
#' n <- 500
#' s <- runif(n, 0.5, 2)
#' x <- rnorm(n, rnorm(n), s)
#' fit_prior <- function(x) {            # captures s, not x
#'   mu  <- mean(x)
#'   tau <- sqrt(max(var(x) - mean(s^2), 0))
#'   function(n) rnorm(n, mu, tau)
#' }
#' lik <- function(theta) rnorm(length(theta), theta, s)   # i-th theta, i-th s
#' eb_calibration_check(x, fit_prior, lik)
#' }
#' @export
eb_calibration_check <- function(X_obs,
                                 fit_prior,
                                 likelihood_sampler,
                                 feature_transform = function(x) x,
                                 B = 199,
                                 K = 5,
                                 alpha = 0.05,
                                 standardize = TRUE,
                                 cores = 1L,
                                 seed = 1,
                                 verbose = TRUE) {

  # ---- Input validation -----------------------------------------------------
  if (missing(fit_prior))
    stop("`fit_prior` is required: the test refits G_hat on every reference ",
         "dataset, so it needs the g-modeling routine, function(x) returning ",
         "a prior sampler, not a sampler from one fitted G_hat.", call. = FALSE)
  if (!is.function(fit_prior))
    stop("`fit_prior` must be a function(x) returning a prior sampler.")
  if (!is.function(likelihood_sampler)) stop("`likelihood_sampler` must be a function.")
  if (!is.function(feature_transform))  stop("`feature_transform` must be a function.")
  B <- check_count(B, "B")
  K <- check_count(K, "K")
  cores <- check_count(cores, "cores")
  check_alpha(alpha)
  if (cores > 1L && .Platform$OS.type == "windows") {
    warning("`cores > 1` is not supported on Windows; running serially.",
            call. = FALSE)
    cores <- 1L
  }
  if (!requireNamespace("energy", quietly = TRUE))
    stop("Package 'energy' is required. Install it with install.packages('energy').")

  X_obs <- as_dataset(X_obs)
  n     <- n_units(X_obs)

  # ---- Protect the caller's random number stream -----------------------------
  # With a seed, the caller's stream is put back exactly as it was. Without
  # one, the test draws from the caller's stream and leaves it advanced past
  # those draws, but the per-replicate reseeding below is still undone.
  rng_on_exit <- get_rng_state()
  on.exit(set_rng_state(rng_on_exit), add = TRUE)
  if (!is.null(seed)) set.seed(seed)

  # ---- Fitted model objects ---------------------------------------------------
  as_model <- function(fit) {
    if (is.function(fit)) {
      list(prior_sampler = fit, likelihood_sampler = likelihood_sampler)
    } else if (is.list(fit) && is.function(fit$prior_sampler)) {
      list(prior_sampler      = fit$prior_sampler,
           likelihood_sampler = if (is.function(fit$likelihood_sampler))
             fit$likelihood_sampler else likelihood_sampler)
    } else {
      stop("`fit_prior` must return a function(n) or a list with a ",
           "`prior_sampler` function.")
    }
  }
  fit_model <- function(x) as_model(fit_prior(x))

  simulate <- function(model) {
    sim <- as_dataset(model$likelihood_sampler(model$prior_sampler(n)))
    if (n_units(sim) != n)
      stop("Simulated dataset has ", n_units(sim), " units; expected ", n, ".")
    sim
  }

  # ---- Feature map --------------------------------------------------------------
  # Scaling constants come from the dataset being tested and are reused for
  # its pseudo-samples, so the same map is applied to every dataset.
  features <- function(x, scaling = NULL) {
    Z <- as.matrix(feature_transform(x))
    storage.mode(Z) <- "double"
    if (!standardize) return(list(Z = Z, scaling = NULL))
    if (is.null(scaling)) {
      mu <- colMeans(Z)
      s  <- apply(Z, 2L, stats::sd)
      s[!is.finite(s) | s == 0] <- 1
      scaling <- list(mu = mu, s = s)
    }
    Z <- sweep(sweep(Z, 2L, scaling$mu, "-"), 2L, scaling$s, "/")
    list(Z = Z, scaling = scaling)
  }

  # ---- The map, applied identically to every dataset (lines 2-5, 9-12) ---------
  phi <- function(x, model) {
    fx <- features(x)
    vapply(seq_len(K), function(k) {
      energy_distance(fx$Z, features(simulate(model), fx$scaling)$Z)
    }, numeric(1))
  }

  # ---- Observed side (lines 1-5) ------------------------------------------------
  model_obs <- fit_model(X_obs)
  T_obs_k   <- phi(X_obs, model_obs)
  T_obs     <- stats::median(T_obs_k)

  # ---- Reference side (lines 6-13) ----------------------------------------------
  # Per-replicate seeds make the result independent of `cores`.
  rep_seeds <- sample.int(.Machine$integer.max, B)
  if (is.null(seed)) rng_on_exit <- get_rng_state()

  one_replicate <- function(b) {
    set.seed(rep_seeds[b])
    tryCatch({
      x_b     <- simulate(model_obs)                # line 7
      model_b <- fit_model(x_b)                     # line 8
      stats::median(phi(x_b, model_b))              # lines 9-12
    }, error = function(e) structure(NA_real_, msg = conditionMessage(e)))
  }

  if (verbose)
    cat(sprintf("EB calibration check: B = %d, K = %d, n = %d\n", B, K, n))

  if (cores > 1L) {
    res <- parallel::mclapply(seq_len(B), one_replicate, mc.cores = cores,
                              mc.set.seed = FALSE)
  } else {
    res <- vector("list", B)
    if (verbose) pb <- utils::txtProgressBar(min = 0L, max = B, style = 3L)
    for (b in seq_len(B)) {
      res[[b]] <- one_replicate(b)
      if (verbose) utils::setTxtProgressBar(pb, b)
    }
    if (verbose) { close(pb); cat("\n") }
  }

  T_ref  <- vapply(res, function(r) if (is.numeric(r)) r[1] else NA_real_, numeric(1))
  failed <- !is.finite(T_ref)
  n_failed <- sum(failed)
  if (n_failed > 0L) {
    msgs <- unique(unlist(lapply(res[failed], attr, "msg")))
    warning(sprintf("%d of %d reference replicates failed. ", n_failed, B),
            "They are counted as exceeding T_obs in `p_value` (conservative). ",
            "First error: ", if (length(msgs)) msgs[1] else "non-finite statistic",
            call. = FALSE)
  }
  T_ref <- T_ref[!failed]
  B_eff <- length(T_ref)
  if (B_eff == 0L) stop("All reference replicates failed.")

  # ---- Monte Carlo p-value (line 14) --------------------------------------------
  # Failures are treated as exceeding T_obs: if the refit fails more often on
  # extreme synthetic datasets, dropping them would remove large T_ref values
  # and make the test anti-conservative.
  n_exceed        <- sum(T_ref >= T_obs)
  p_value         <- (1 + n_exceed + n_failed) / (B + 1)
  p_value_dropped <- (1 + n_exceed) / (B_eff + 1)

  structure(
    list(
      p_value         = p_value,
      p_value_dropped = p_value_dropped,
      T_obs           = T_obs,
      T_obs_k         = T_obs_k,
      T_ref           = T_ref,
      B               = B,
      B_eff           = B_eff,
      K               = K,
      n_failed        = n_failed,
      alpha           = alpha,
      n               = n
    ),
    class = "eb_calibration"
  )
}


# ---- RNG state helpers ----------------------------------------------------------

# NULL when the caller has not yet initialised the generator.
get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

set_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
}

check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  invisible(alpha)
}

# ---- Internal helpers ---------------------------------------------------------

check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != round(x))
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  as.integer(x)
}

# Canonical dataset form: numeric vector for one column, numeric matrix otherwise.
as_dataset <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x) && ncol(x) == 1L) x <- x[, 1L]
  if (!is.matrix(x)) x <- as.vector(x)
  x
}

n_units <- function(x) if (is.matrix(x)) nrow(x) else length(x)

# Two-sample energy distance, V-statistic form (proposal eq. 6).
energy_distance <- function(X, Y) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (ncol(X) == 1L) return(energy_distance_1d(X[, 1L], Y[, 1L]))
  nx <- nrow(X); ny <- nrow(Y)
  e  <- energy::eqdist.e(rbind(X, Y), sizes = c(nx, ny))
  unname(e) * (nx + ny) / (nx * ny)
}

# Exact O(n log n) energy distance for univariate samples, via sorting:
# sum_{i<j} |s_i - s_j| = sum_i (2i - m - 1) s_(i).
energy_distance_1d <- function(x, y) {
  half_sum <- function(s) {
    s <- sort(s); m <- length(s)
    sum((2 * seq_len(m) - m - 1) * s)
  }
  nx <- length(x); ny <- length(y)
  Wx <- half_sum(x); Wy <- half_sum(y)
  cross <- half_sum(c(x, y)) - Wx - Wy        # sum_{i,j} |x_i - y_j|
  2 * cross / (nx * ny) - 2 * Wx / nx^2 - 2 * Wy / ny^2
}



# ---- S3 methods -----------------------------------------------------------------

#' Print an EB calibration check
#'
#' @param x An `"eb_calibration"` object.
#' @param alpha Level at which to report a rejection. Defaults to the level
#'   the test was run with.
#' @param ... Ignored.
#' @export
print.eb_calibration <- function(x, alpha = x$alpha, ...) {
  if (is.null(alpha)) alpha <- 0.05
  check_alpha(alpha)
  cat("Empirical Bayes Goodness-of-Fit Calibration Check\n")
  cat("--------------------------------------------------\n")
  cat(sprintf("  n                    : %d\n", x$n))
  cat(sprintf("  Reference datasets B : %d\n", x$B))
  if (x$n_failed > 0)
    cat(sprintf("  Failed refits        : %d (counted as exceeding T_obs)\n",
                x$n_failed))
  cat(sprintf("  Inner draws (K)      : %d\n", x$K))
  cat(sprintf("  T_obs (median of K)  : %.4g\n", x$T_obs))
  cat(sprintf("  Monte Carlo p-value  : %.4f\n", x$p_value))
  if (x$n_failed > 0)
    cat(sprintf("  p (failures dropped) : %.4f\n", x$p_value_dropped))
  cat("\n")
  if (x$p_value <= alpha) {
    cat(sprintf("  >> Reject at level %g: the fitted prior family does not\n", alpha))
    cat("     reproduce the observed marginal distribution.\n")
  } else {
    cat(sprintf("  >> No rejection at level %g: the fitted prior family is\n", alpha))
    cat("     consistent with the observed marginal. This does not certify\n")
    cat("     the likelihood or independence across units.\n")
  }
  invisible(x)
}


#' @export
plot.eb_calibration <- function(x, ...) {
  rng <- range(c(x$T_ref, x$T_obs, x$T_obs_k))
  graphics::hist(x$T_ref,
                 breaks      = 30L,
                 col         = "steelblue",
                 border      = "white",
                 xlim        = rng,
                 xlab        = "Energy distance (median of K)",
                 main        = "Reference distribution vs T_obs",
                 probability = TRUE,
                 ...)
  graphics::abline(v = x$T_obs, col = "firebrick", lwd = 2L, lty = 2L)
  graphics::rug(x$T_obs_k, col = "firebrick")
  graphics::legend("topright",
                   legend = c(expression(T[b]^{ref}),
                              bquote(T^{obs} == .(signif(x$T_obs, 3))),
                              bquote(p == .(round(x$p_value, 4)))),
                   col = c("steelblue", "firebrick", NA),
                   lwd = c(8L, 2L, NA), lty = c(1L, 2L, NA),
                   bty = "n")
  invisible(x)
}
