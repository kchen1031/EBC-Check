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
#' @section Parallel computing on Windows:
#' Windows cannot fork, so with `cores > 1` each worker is a new R session.
#' Before the refits start, each worker gets the caller's library paths and
#' random number generator kind, loads this package (the development copy,
#' via `pkgload::load_all()`, if the package was loaded with
#' `devtools::load_all()`), attaches `packages`, and receives the objects
#' named in `export`. Closures returned by [get_data()] carry their own data
#' and need nothing exported. Closures written at top level that use other
#' top-level objects do: pass those names in `export`. If the first batch of
#' refits fails because something is missing on the workers, the function
#' stops at once with the error rather than running all `B` refits. Windows
#' Firewall may ask once whether R may open local network connections; the
#' workers only talk to the current session on this machine.
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
#' @param cores Number of worker processes for the `B` refits. Default 1
#'   (serial). On macOS and Linux the workers are forked with
#'   [parallel::mclapply()] and see everything in the current session. On
#'   Windows they are fresh R sessions in a socket cluster
#'   ([parallel::makePSOCKcluster()]); see the "Parallel computing on Windows"
#'   section for what is copied to them. Results are identical for any value
#'   of `cores`, on any platform, given the same `seed`. Values above `B` are
#'   reduced to `B`.
#' @param export Character vector of object names, looked up from where
#'   `eb_calibration_check()` is called, to copy to socket-cluster workers
#'   (Windows only). Needed when `fit_prior`, `likelihood_sampler` or
#'   `feature_transform` were defined at top level and use other top-level
#'   objects, such as a vector of standard errors `s`. Ignored when the
#'   workers are forked. Default `NULL`.
#' @param packages Character vector of packages to attach on socket-cluster
#'   workers (Windows only). Default `NULL` attaches every non-base package
#'   attached in the current session, so unqualified calls such as `GLmix()`
#'   after `library(REBayes)` work on the workers. Use `character(0)` to
#'   attach none. Ignored when the workers are forked.
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
#'
#' # Same test on 4 cores. On Windows, `s` must be exported because both
#' # closures look it up at top level; elsewhere `export` is ignored.
#' eb_calibration_check(x, fit_prior, lik, cores = 4, export = "s")
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
                                 export = NULL,
                                 packages = NULL,
                                 seed = 1,
                                 verbose = TRUE) {

  # Where `export` names are looked up. Removed again before the replicate
  # function is sent to workers, so a large calling frame is not copied.
  export_env <- parent.frame()

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
  cores <- min(check_count(cores, "cores"), B)
  check_alpha(alpha)
  if (!is.null(export) && !is.character(export))
    stop("`export` must be NULL or a character vector of object names.",
         call. = FALSE)
  if (!is.null(packages) && !is.character(packages))
    stop("`packages` must be NULL or a character vector of package names.",
         call. = FALSE)
  backend <- parallel_backend(cores)
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
    cat(sprintf("EB calibration check: B = %d, K = %d, n = %d%s\n", B, K, n,
                if (cores > 1L) sprintf(", cores = %d (%s)", cores, backend)
                else ""))

  # ---- Run the B replicates --------------------------------------------------
  # `run_batch(idx)` returns the replicates for indices `idx`, in order.
  if (backend == "psock") {
    if (verbose) cat("Starting worker processes...\n")
    cl <- start_psock_cluster(cores, packages)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (length(export)) {
      missing_obj <- export[!vapply(export, exists, logical(1),
                                    envir = export_env)]
      if (length(missing_obj))
        stop("`export` names objects that do not exist: ",
             paste(missing_obj, collapse = ", "), call. = FALSE)
      parallel::clusterExport(cl, export, envir = export_env)
    }
    rm(export_env)
    parallel::clusterCall(cl, worker_set, one_replicate)
    run_batch <- function(idx) parallel::parLapplyLB(cl, idx, worker_run)
  } else if (backend == "fork") {
    run_batch <- function(idx)
      parallel::mclapply(idx, one_replicate, mc.cores = cores,
                         mc.set.seed = FALSE)
  } else {
    run_batch <- function(idx) lapply(idx, one_replicate)
  }

  # Serial: one replicate per batch. Parallel: a first batch of `cores`
  # replicates, which catches set-up problems early; then, if printing
  # progress, batches of whole multiples of `cores` (about 20 updates),
  # otherwise everything else at once for the best load balancing.
  batches <- if (backend == "serial") {
    as.list(seq_len(B))
  } else {
    rest <- if (B > cores) seq.int(cores + 1L, B) else integer(0)
    size <- if (verbose) cores * max(1L, round(B / (20 * cores))) else length(rest)
    c(list(seq_len(cores)),
      if (length(rest)) split(rest, ceiling(seq_along(rest) / size)))
  }

  res <- vector("list", B)
  if (verbose) pb <- utils::txtProgressBar(min = 0L, max = B, style = 3L)
  for (i in seq_along(batches)) {
    idx <- batches[[i]]
    res[idx] <- run_batch(idx)
    if (i == 1L && backend == "psock") check_worker_setup(res[idx])
    if (verbose) utils::setTxtProgressBar(pb, max(idx))
  }
  if (verbose) { close(pb); cat("\n") }

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

# ---- Parallel helpers -----------------------------------------------------------

# "serial", "fork" (macOS / Linux) or "psock" (Windows). The internal option
# `ebcalibration.backend = "psock"` forces socket workers on any platform, so
# the Windows path can be tested elsewhere.
parallel_backend <- function(cores) {
  if (cores == 1L) return("serial")
  if (.Platform$OS.type == "windows" ||
      identical(getOption("ebcalibration.backend"), "psock")) "psock" else "fork"
}

# Start `cores` socket workers set up like the current session. Functions sent
# to workers here are base or pkgload functions, never closures from this
# package, because this package's namespace is not loaded on a worker until
# the last step.
start_psock_cluster <- function(cores, packages = NULL) {
  cl <- parallel::makePSOCKcluster(cores)
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl), add = TRUE)

  # `.libPaths` keeps its paths in its own enclosing environment, so sending
  # the function itself would only change a copy. Look it up on the worker.
  set_lib_paths <- function(p) { .libPaths(p); NULL }
  environment(set_lib_paths) <- baseenv()
  parallel::clusterCall(cl, set_lib_paths, .libPaths())
  kind <- RNGkind()
  parallel::clusterCall(cl, RNGkind, kind[1], kind[2], kind[3])

  # This package: the development copy if loaded with devtools::load_all(),
  # otherwise the installed copy, which loads when the first closure arrives.
  ns <- environment(eb_calibration_check)
  if (exists(".__DEVTOOLS__", envir = ns, inherits = FALSE)) {
    if (!requireNamespace("pkgload", quietly = TRUE))
      stop("The package was loaded with load_all() but 'pkgload' is not ",
           "installed, so workers cannot load it.", call. = FALSE)
    parallel::clusterCall(cl, getExportedValue("pkgload", "load_all"),
                          getNamespaceInfo(ns, "path"),
                          export_all = FALSE, helpers = FALSE, quiet = TRUE)
  } else {
    parallel::clusterCall(cl, loadNamespace, getNamespaceName(ns))
  }

  # Attach packages in their current search-path order (last attached ends
  # up first), skipping base packages and this package.
  if (is.null(packages)) {
    base_pkgs <- c("base", "methods", "datasets", "utils", "grDevices",
                   "graphics", "stats")
    packages <- rev(setdiff(.packages(), c(base_pkgs, getNamespaceName(ns))))
  }
  failed <- character(0)
  for (p in packages) {
    attached <- tryCatch({
      parallel::clusterCall(cl, library, p, character.only = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (!attached) failed <- c(failed, p)
  }
  if (length(failed))
    warning("Could not attach on workers: ", paste(failed, collapse = ", "),
            call. = FALSE)

  ok <- TRUE
  cl
}

# Worker-side storage for the replicate function, so it is sent to each worker
# once rather than with every task. These two functions live in this package's
# namespace, which is serialized by reference and resolved on the worker.
.worker_state <- new.env(parent = emptyenv())
worker_set <- function(f) {
  assign("replicate", f, envir = .worker_state)
  NULL
}
worker_run <- function(b) get("replicate", envir = .worker_state)(b)

# Stop early if every replicate in the first parallel batch failed because
# something did not reach the workers.
check_worker_setup <- function(res) {
  msgs <- vapply(res, function(r) {
    m <- attr(r, "msg")
    if (is.null(m)) NA_character_ else m
  }, character(1))
  setup_error <- grepl("not found|could not find function|there is no package",
                       msgs)
  if (length(msgs) && all(setup_error)) {
    stop("Every refit in the first parallel batch failed on the workers:\n  ",
         msgs[1], "\n",
         "Workers on Windows are new R sessions. If `fit_prior`, ",
         "`likelihood_sampler` or `feature_transform` use top-level objects, ",
         "pass their names in `export`; if they call functions from a package ",
         "without `pkg::`, list it in `packages`.", call. = FALSE)
  }
  invisible(NULL)
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
