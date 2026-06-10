#' Get dataset, fitted G_hat, sampler functions, and feature transform
#'
#' @param paper Character string identifying the paper (see Details)
#' @param dataset Character string identifying the dataset within the paper
#'
#' @return A list with elements:
#' \describe{
#'   \item{data}{The observed dataset as a data frame}
#'   \item{X_obs}{The observed statistic(s) to pass to eb_calibration_check,
#'         on the same scale as likelihood_sampler output (vector or matrix)}
#'   \item{G_hat}{The fitted empirical Bayes prior (structure varies by model)}
#'   \item{prior_sampler}{function(n): draws n samples of theta from G_hat}
#'   \item{likelihood_sampler}{function(theta): draws one pseudo-observation
#'         per element of theta, on the NATURAL scale (raw counts / values)}
#'   \item{feature_transform}{function(x): maps natural-scale observations to
#'         the scale on which the energy statistic is computed}
#' }
#'
#' @examples
#' \dontrun{
#' bundle <- get_data("Willwerscheid", "wOBA")
#' result <- eb_calibration_check(
#'   X_obs              = bundle$X_obs,
#'   prior_sampler      = bundle$prior_sampler,
#'   likelihood_sampler = bundle$likelihood_sampler,
#'   feature_transform  = bundle$feature_transform
#' )
#' }
#'
#' @export
get_data <- function(paper, dataset) {

  key <- paste0(tolower(trimws(paper)), "_", tolower(trimws(dataset)))

  # ------------------------------------------------------------------
  # Willwerscheid et al. 2024 - wOBA (ebnm package)
  #   x_i | theta_i ~ N(theta_i, s_i^2),  theta_i ~ G
  #   G fit by one of three families (paper Section 5):
  #     "normal"   -> ebnm(..., prior_family = "normal",   mode = "estimate")
  #     "unimodal" -> ebnm(..., prior_family = "unimodal", mode = "estimate")
  #     "npmle"    -> ebnm(..., prior_family = "npmle")
  # Keys: willwerscheid_woba_normal, willwerscheid_woba_unimodal,
  #       willwerscheid_woba_npmle
  # ------------------------------------------------------------------
  if (key %in% c("willwerscheid_woba_normal",
                 "willwerscheid_woba_unimodal",
                 "willwerscheid_woba_npmle")) {

    if (!requireNamespace("ebnm", quietly = TRUE))
      stop("Please install the 'ebnm' package.")

    data("wOBA", package = "ebnm", envir = environment())
    s_obs <- wOBA$s

    family <- sub("^willwerscheid_woba_", "", key)   # normal / unimodal / npmle

    fit <- if (family == "npmle") {
      ebnm::ebnm(x = wOBA$x, s = wOBA$s, prior_family = "npmle")
    } else {
      ebnm::ebnm(x = wOBA$x, s = wOBA$s,
                 prior_family = family, mode = "estimate")
    }
    G_hat <- fit$fitted_g

    # Generic draw from a fitted ebnm prior:
    #   normalmix (normal, npmle) -> $mean/$sd  (npmle has sd == 0)
    #   unimix    (unimodal)      -> $a/$b      (uniform mixture)
    sample_from_g <- function(g, n) {
      idx <- sample.int(length(g$pi), size = n, replace = TRUE, prob = g$pi)
      if (!is.null(g$mean) && !is.null(g$sd)) {
        stats::rnorm(n, mean = g$mean[idx], sd = g$sd[idx])
      } else if (!is.null(g$a) && !is.null(g$b)) {
        stats::runif(n, min = g$a[idx], max = g$b[idx])
      } else {
        stop("Unsupported fitted_g structure for wOBA prior sampler.")
      }
    }

    return(list(
      data  = wOBA,
      X_obs = wOBA$x,
      G_hat = G_hat,

      prior_sampler = function(n) sample_from_g(G_hat, n),

      likelihood_sampler = function(theta) {
        n <- length(theta)
        s <- rep_len(s_obs, n)
        stats::rnorm(n, mean = theta, sd = s)
      },

      feature_transform = function(x) x
    ))
  }

  # ------------------------------------------------------------------
  # 2. Koenker et al. 2016 - Norberg (REBayes package)
  # ------------------------------------------------------------------
  if (key == "koenker_norberg") {

    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install the 'REBayes' package.")

    data("Norberg", package = "REBayes", envir = environment())
    x_obs <- Norberg$Death
    e_obs <- Norberg$Exposure / 344

    fit   <- REBayes::Pmix(x = x_obs, v = e_obs)
    G_hat <- list(support = fit$x, weights = fit$y)

    return(list(
      data  = Norberg,
      X_obs = x_obs,
      G_hat = G_hat,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights)
      },

      likelihood_sampler = function(theta) {
        n      <- length(theta)
        e      <- rep_len(e_obs, n)
        stats::rpois(n, lambda = theta * e)
      },

      # Bivariate feature: (log exposure, log1p count). Returns a 2-column
      # matrix so the energy test compares the JOINT (E, X) distribution.
      feature_transform = function(x) {
        x <- as.numeric(x)
        e <- rep_len(e_obs, length(x))
        cbind(
          logE = log(pmax(e, 1e-12)),
          rate = log1p(x)
        )
      }
    ))
  }

  # ------------------------------------------------------------------
  # 3. Koenker et al. 2016 - Flies (REBayes package)
  # ------------------------------------------------------------------
  if (key == "koenker_flies") {

    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install the 'REBayes' package.")

    data("flies", package = "REBayes", envir = environment())
    t_obs <- flies$age

    fit   <- REBayes::WGLVmix(t_obs)
    G_hat <- list(support = fit$x, weights = fit$y, alpha = fit$alpha)

    return(list(
      data  = flies,
      X_obs = t_obs,
      G_hat = G_hat,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights)
      },

      likelihood_sampler = function(theta) {
        alpha <- G_hat$alpha
        scale <- theta^(-1 / alpha)
        stats::rweibull(length(theta), shape = alpha, scale = scale)
      },

      feature_transform = function(x) x
    ))
  }

  # ------------------------------------------------------------------
  # 4. Narasimhan et al. 2020 - Surgery (deconvolveR package)
  # ------------------------------------------------------------------
  if (key == "narasimhan_surgery") {

    if (!requireNamespace("deconvolveR", quietly = TRUE))
      stop("Please install the 'deconvolveR' package.")

    data("surg", package = "deconvolveR", envir = environment())
    n_obs <- surg$n
    x_obs <- surg$s

    tau <- seq(0.01, 0.99, length.out = 200)

    X_mat <- cbind(surg$n, surg$s)
    fit   <- deconvolveR::deconv(tau = tau, X = X_mat, family = "Binomial")

    G_hat <- list(support = tau, weights = fit$stats[, "g"])

    return(list(
      data  = surg,
      X_obs = x_obs,
      G_hat = G_hat,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },

      likelihood_sampler = function(theta) {
        n <- length(theta)
        stats::rbinom(n, size = rep_len(n_obs, n), prob = theta)
      },

      feature_transform = function(x) {
        n <- rep_len(n_obs, length(x))
        asin(sqrt(x / pmax(n, 1e-12)))
      }
    ))
  }

  # ------------------------------------------------------------------
  # 5. Narasimhan et al. 2020 - Shakespeare (deconvolveR package)
  # ------------------------------------------------------------------
  if (key == "narasimhan_shakespeare") {

    if (!requireNamespace("deconvolveR", quietly = TRUE))
      stop("Please install the 'deconvolveR' package.")

    data("bardWordCount", package = "deconvolveR", envir = environment())

    x_obs <- rep(seq_along(bardWordCount), times = bardWordCount)

    lambda <- seq(-4, 4.5, 0.025)
    tau    <- exp(lambda)
    fit    <- deconvolveR::deconv(tau = tau, y = bardWordCount,
                                  n = 100, family = "Poisson", c0 = 2)
    G_hat  <- list(support = tau, weights = fit$stats[, "g"])

    return(list(
      data  = bardWordCount,
      X_obs = x_obs,
      G_hat = G_hat,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },

      likelihood_sampler = function(theta) {
        n   <- length(theta)
        out <- integer(n)
        for (i in seq_len(n)) {
          repeat {
            x <- stats::rpois(1, lambda = theta[i])
            if (x >= 1) { out[i] <- x; break }
          }
        }
        out
      },

      feature_transform = function(x) sqrt(x)
    ))
  }

  # ------------------------------------------------------------------
  # Gu & Koenker (Bayesball) - location-scale mixture (REBayes::bball)
  # ------------------------------------------------------------------
  if (key == "gu_baseball") {

    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install 'REBayes' (and MOSEK via Rmosek).")

    data("bball", package = "REBayes", envir = environment())

    # --- Section 4.1 filtering ---
    bb <- bball
    bb <- bb[bb$AB > 10, , drop = FALSE]
    bb <- bb[bb$pitcher == 0, , drop = FALSE]
    hs <- table(bb$id)
    bb <- bb[bb$id %in% names(hs)[hs >= 3], , drop = FALSE]

    # --- Transform + per-player sufficient statistics ---
    y   <- asin(sqrt((bb$H + 0.25) / (bb$AB + 0.5)))
    nu2 <- 1 / (4 * bb$AB)

    players <- unique(bb$id)
    np      <- length(players)
    muhat <- S <- v2 <- r <- numeric(np)
    for (j in seq_len(np)) {
      sel  <- bb$id == players[j]
      yj   <- y[sel]; nu2j <- nu2[sel]; mj <- length(yj)
      w    <- 1 / nu2j
      muhat[j] <- sum(w * yj) / sum(w)
      S[j]     <- sum((yj - muhat[j])^2 / nu2j) / (mj - 1)
      v2[j]    <- 1 / (4 * sum(bb$AB[sel]))
      r[j]     <- (mj - 1) / 2
    }

    # --- Fit the JOINT (mu, theta) NPMLE ---
    fit <- REBayes::WGLVmix(y = y, id = bb$id, w = 1 / nu2, u = 100, v = 120)
    grid_mu    <- fit$u
    grid_theta <- fit$v
    W          <- fit$fuv

    G <- expand.grid(mu = grid_mu, theta = grid_theta)
    G$w <- as.vector(W)
    G   <- G[G$w > 1e-10, ]
    G$w <- G$w / sum(G$w)

    baseball_ls <- data.frame(id = players, muhat = muhat, S = S)

    return(list(
      data  = baseball_ls,
      X_obs = cbind(muhat = muhat, S = S),
      G_hat = G,

      prior_sampler = function(n) {
        idx <- sample(nrow(G), size = n, replace = TRUE, prob = G$w)
        cbind(mu = G$mu[idx], theta = G$theta[idx])
      },

      likelihood_sampler = function(theta_mat) {
        mu <- theta_mat[, 1]; th <- theta_mat[, 2]
        cbind(
          muhat = stats::rnorm(length(mu), mean = mu, sd = sqrt(th * v2)),
          S     = stats::rgamma(length(th), shape = r, scale = th / r)
        )
      },

      feature_transform = function(x) x
    ))
  }

  # ============================================================
  # Species-richness datasets (Baek & Park 2022)
  #   x_i | lambda_i ~ Poisson(lambda_i), observed only if x_i >= 1
  #   (zero-truncated), lambda_i ~ G.
  # Keys: baek_brainvessel, baek_hardcandy, baek_deathnotice,
  #       baek_traffic, baek_microbial, baek_shakespeare
  # ============================================================
  if (startsWith(key, "baek_")) {

    fit_ztpois_mixture <- function(x, n_grid = 100, max_iter = 2000, tol = 1e-9) {
      x <- as.numeric(x)
      lam_max <- max(2 * max(x), 5)
      grid    <- exp(seq(log(0.05), log(lam_max), length.out = n_grid))
      K       <- length(grid)

      logA <- outer(x, grid, function(xi, lam) {
        stats::dpois(xi, lam, log = TRUE) - log1p(-exp(-lam))
      })

      w <- rep(1 / K, K)
      for (iter in seq_len(max_iter)) {
        logpost <- sweep(logA, 2, log(w + 1e-300), "+")
        mx <- apply(logpost, 1, max)
        P  <- exp(logpost - mx)
        P  <- P / rowSums(P)
        w_new <- colMeans(P)
        w_new <- w_new / sum(w_new)
        if (max(abs(w_new - w)) < tol) { w <- w_new; break }
        w <- w_new
      }
      keep <- w > 1e-6
      list(grid = grid[keep], w = w[keep] / sum(w[keep]))
    }

    species_registry <- list(
      brainvessel = list(counts = 1:14,
                         freqs = c(4,15,31,39,55,54,49,47,31,16,9,8,4,3), true_c = 366),
      hardcandy   = list(counts = 1:20,
                         freqs = c(54,49,62,44,25,26,15,15,10,10,10,10,3,3,5,5,4,1,2,1), true_c = 456),
      deathnotice = list(counts = 1:9,
                         freqs = c(267,271,185,111,61,27,8,3,1), true_c = 1096),
      traffic     = list(counts = 1:7,
                         freqs = c(1317,239,42,14,4,4,1), true_c = 9461),
      microbial   = list(counts = c(1,2,3,4,5,6,7,9,11,13,14,16,21,27,32,43),
                         freqs = c(381,65,23,18,4,5,3,1,4,3,2,1,1,1,1,1), true_c = NA)
    )

    ds <- sub("^baek_", "", key)

    if (ds == "shakespeare") {
      if (!requireNamespace("deconvolveR", quietly = TRUE))
        stop("Please install 'deconvolveR' for the Shakespeare data.")
      data("bardWordCount", package = "deconvolveR", envir = environment())
      reg <- list(counts = seq_along(bardWordCount),
                  freqs  = as.numeric(bardWordCount), true_c = NA)
    } else {
      if (!ds %in% names(species_registry))
        stop("Unknown species dataset: '", ds, "'. Options: ",
             paste(names(species_registry), collapse = ", "), ", shakespeare")
      reg <- species_registry[[ds]]
    }

    x_obs <- rep(reg$counts, times = reg$freqs)
    stopifnot(min(x_obs) >= 1)

    fit   <- fit_ztpois_mixture(x_obs, n_grid = 100)
    G_hat <- list(support = fit$grid, weights = fit$w)

    return(list(
      data   = data.frame(count = reg$counts, freq = reg$freqs),
      X_obs  = x_obs,
      G_hat  = G_hat,
      true_c = reg$true_c,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },

      likelihood_sampler = function(theta) {
        n  <- length(theta)
        p0 <- stats::dpois(0, theta)
        stats::qpois(stats::runif(n, p0, 1), theta)
      },

      feature_transform = function(x) sqrt(x)
    ))
  }


  # ------------------------------------------------------------------
  # Jana et al. (mindist) - Hockey (NHL goals, Poisson EB)
  # Keys: jana_hockey_npmle, jana_hockey_hellinger, jana_hockey_chisquare
  #
  # Model:  Y_i | theta_i ~ Poisson(theta_i),  theta_i ~ G
  # G is fit from 2017-18 goals via one of three minimum-distance
  # estimators (NPMLE / Hellinger / chi-square). 2018-19 goals are
  # stored in the data slot for downstream prediction checks.
  #
  # Data: hockey-reference.com, skaters appearing in BOTH seasons
  #       (matched on player name), per Jana et al. (2024).
  # ------------------------------------------------------------------
  if (key %in% c("jana_hockey_npmle", "jana_hockey_hellinger", "jana_hockey_chisquare")) {

    if (!requireNamespace("rvest", quietly = TRUE))
      stop("Please install the 'rvest' package.")
    if (!requireNamespace("dplyr", quietly = TRUE))
      stop("Please install the 'dplyr' package.")

    `%>%` <- dplyr::`%>%`

    # --- 1. scrape, or reuse a pre-built data frame via options() ------
    dat_opt <- getOption("ebcalibration.hockey_dat")
    if (!is.null(dat_opt)) {
      dat <- dat_opt
    } else {
      url18 <- "https://www.hockey-reference.com/leagues/NHL_2018_skaters.html"
      url19 <- "https://www.hockey-reference.com/leagues/NHL_2019_skaters.html"

      sk18 <- rvest::html_table(rvest::read_html(url18), fill = TRUE)[[1]]
      sk19 <- rvest::html_table(rvest::read_html(url19), fill = TRUE)[[1]]

      # Column layout: 1 Rk | 2 Player | 3 Age | 4 Tm | 5 Pos | 6 GP | 7 G | ...
      clean_hockey <- function(df, goal_name) {
        df %>%
          dplyr::select(Player = 2, Team = 4, GP = 6, G = 7) %>%
          dplyr::filter(Player != "Player", !is.na(Player), Player != "") %>%
          dplyr::mutate(
            GP = suppressWarnings(as.numeric(GP)),
            G  = suppressWarnings(as.numeric(G))
          ) %>%
          dplyr::filter(!is.na(GP), !is.na(G)) %>%
          dplyr::group_by(Player) %>%
          dplyr::filter(if (any(Team == "TOT")) Team == "TOT" else TRUE) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup() %>%
          dplyr::select(Player, G) %>%
          dplyr::rename(!!goal_name := G)
      }

      d18 <- clean_hockey(sk18, "goals_2017_18")
      d19 <- clean_hockey(sk19, "goals_2018_19")
      dat <- dplyr::inner_join(d18, d19, by = "Player")
    }

    common <- dat$Player
    y18    <- dat$goals_2017_18   # fit G on these
    y19    <- dat$goals_2018_19   # stored for prediction checks

    # --- 2. minimum-distance Poisson mixture fitter --------------------
    # Minimises sum_y div(p_emp(y), f_G(y)) over discrete G on a grid.
    #   "npmle"     -> KL / standard EM
    #   "hellinger" -> squared Hellinger / multiplicative EM
    #   "chisquare" -> chi^2 / multiplicative EM
    fit_poisson_mindist <- function(y, n_grid = 300, method = "npmle",
                                    max_iter = 5000, tol = 1e-10) {
      y_vals <- sort(unique(y))
      p_emp  <- tabulate(match(y, y_vals)) / length(y)
      h      <- max(2 * max(y_vals), 10)
      grid   <- seq(h / n_grid, h, length.out = n_grid)

      A <- outer(y_vals, grid, stats::dpois)   # nY x K
      w <- rep(1 / n_grid, n_grid)

      for (iter in seq_len(max_iter)) {
        f_mix <- pmax(as.vector(A %*% w), 1e-300)

        w_new <- switch(method,
                        npmle     = w * as.vector(t(A) %*% (p_emp / f_mix)),
                        hellinger = w * as.vector(t(A) %*% sqrt(p_emp / f_mix)),
                        chisquare = w * sqrt(as.vector(t(A) %*% (p_emp^2 / f_mix^2)))
        )

        w_new <- pmax(w_new, 0)
        s     <- sum(w_new)
        if (s < 1e-300) break
        w_new <- w_new / s
        if (max(abs(w_new - w)) < tol) { w <- w_new; break }
        w <- w_new
      }

      keep <- w > 1e-8
      list(support = grid[keep], weights = w[keep] / sum(w[keep]))
    }

    method_str <- sub("^jana_hockey_", "", key)   # npmle / hellinger / chisquare
    fit   <- fit_poisson_mindist(y18, n_grid = 300, method = method_str)
    G_hat <- list(support = fit$support, weights = fit$weights)

    return(list(
      data  = data.frame(player = common, goals18 = y18, goals19 = y19),
      X_obs = y18,
      G_hat = G_hat,

      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },

      likelihood_sampler = function(theta) {
        stats::rpois(length(theta), lambda = theta)
      },

      feature_transform = function(x) sqrt(x)
    ))
  }

  # ------------------------------------------------------------------
  # You (2026) - Project STAR  (Gaussian EB, independent-subgroup)
  # ------------------------------------------------------------------
  # Each STAR school i is one prior "study"; for each of 4 strata g
  # (race x free-lunch cross) we form the within-school small-vs-(reg+aide)
  # difference-in-means estimate theta_hat_ig with known variance sigma2_ig:
  #     theta_hat_ig | theta_ig ~ N(theta_ig, sigma2_ig),  theta_ig ~ G_g.
  # Gaussian EB: G_g = N(tau_g, V_g) fit by profile marginal likelihood.
  # We pool all OBSERVED (i, g) estimates for the marginal energy check.
  # ------------------------------------------------------------------
  if (key == "you_star_gaussian") {

    if (!requireNamespace("AER", quietly = TRUE))
      stop("Please install the 'AER' package.")

    data("STAR", package = "AER", envir = environment())
    G <- 4L

    star <- STAR
    keep <- !is.na(star$readk) & !is.na(star$stark) &
      !is.na(star$lunchk) & !is.na(star$ethnicity) & !is.na(star$schoolidk)
    star <- star[keep, ]

    nonwhite <- as.character(star$ethnicity) != "cauc"
    freelnch <- as.character(star$lunchk) == "free"
    star$stratum <- ifelse( nonwhite &  freelnch, 1L,
                            ifelse( nonwhite & !freelnch, 2L,
                                    ifelse(!nonwhite &  freelnch, 3L, 4L)))
    star$D <- as.integer(as.character(star$stark) == "small")  # small=1, reg/reg+aide=0

    # School x stratum difference-in-means estimates (>=2 obs per arm)
    schools <- unique(star$schoolidk)
    rows <- list()
    for (sch in schools) {
      sub <- star[star$schoolidk == sch, ]
      for (g in seq_len(G)) {
        sg <- sub[sub$stratum == g, ]
        y1 <- sg$readk[sg$D == 1]; y0 <- sg$readk[sg$D == 0]
        if (length(y1) < 2 || length(y0) < 2) next
        rows[[length(rows) + 1L]] <- data.frame(
          school = sch, stratum = g,
          theta_hat = mean(y1) - mean(y0),
          sigma2    = stats::var(y1) / length(y1) + stats::var(y0) / length(y0)
        )
      }
    }
    est <- do.call(rbind, rows)

    # Per-stratum Gaussian EB via profile marginal likelihood
    #   theta_hat ~ N(tau_g, sigma2 + V_g)
    tau_g <- numeric(G); V_g <- numeric(G)
    for (g in seq_len(G)) {
      eg  <- est[est$stratum == g, ]
      yg  <- eg$theta_hat; s2g <- eg$sigma2
      nll <- function(par) {
        tv <- s2g + exp(par[2])
        0.5 * sum(log(tv) + (yg - par[1])^2 / tv)
      }
      init <- c(mean(yg), log(max(stats::var(yg) - mean(s2g), 1)))
      fit  <- stats::optim(init, nll, method = "Nelder-Mead",
                           control = list(maxit = 2000, reltol = 1e-10))
      tau_g[g] <- fit$par[1]; V_g[g] <- exp(fit$par[2])
    }

    strata_obs <- est$stratum
    s2_obs     <- est$sigma2
    counts_g   <- as.numeric(table(factor(strata_obs, levels = seq_len(G))))

    return(list(
      data  = est,
      X_obs = est$theta_hat,
      G_hat = list(tau = tau_g, V = V_g, strata = seq_len(G)),

      prior_sampler = function(n) {
        g <- sample(seq_len(G), size = n, replace = TRUE, prob = counts_g)
        stats::rnorm(n, mean = tau_g[g], sd = sqrt(V_g[g]))
      },

      likelihood_sampler = function(theta) {
        n  <- length(theta)
        s2 <- sample(s2_obs, size = n, replace = TRUE)
        stats::rnorm(n, mean = theta, sd = sqrt(s2))
      },

      feature_transform = function(x) x
    ))
  }

  # ------------------------------------------------------------------
  # You (2026) - Project STAR  (NPMLE EB, independent-subgroup)
  # ------------------------------------------------------------------
  # Same model as you_star_gaussian, but G_g is the Kiefer-Wolfowitz NPMLE
  # (Gaussian location mixture) fit per stratum via REBayes::GLmix.
  # ------------------------------------------------------------------
  if (key == "you_star_npmle") {

    if (!requireNamespace("AER", quietly = TRUE))
      stop("Please install the 'AER' package.")
    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install the 'REBayes' package.")

    data("STAR", package = "AER", envir = environment())
    G <- 4L

    star <- STAR
    keep <- !is.na(star$readk) & !is.na(star$stark) &
      !is.na(star$lunchk) & !is.na(star$ethnicity) & !is.na(star$schoolidk)
    star <- star[keep, ]

    nonwhite <- as.character(star$ethnicity) != "cauc"
    freelnch <- as.character(star$lunchk) == "free"
    star$stratum <- ifelse( nonwhite &  freelnch, 1L,
                            ifelse( nonwhite & !freelnch, 2L,
                                    ifelse(!nonwhite &  freelnch, 3L, 4L)))
    star$D <- as.integer(as.character(star$stark) == "small")

    schools <- unique(star$schoolidk)
    rows <- list()
    for (sch in schools) {
      sub <- star[star$schoolidk == sch, ]
      for (g in seq_len(G)) {
        sg <- sub[sub$stratum == g, ]
        y1 <- sg$readk[sg$D == 1]; y0 <- sg$readk[sg$D == 0]
        if (length(y1) < 2 || length(y0) < 2) next
        rows[[length(rows) + 1L]] <- data.frame(
          school = sch, stratum = g,
          theta_hat = mean(y1) - mean(y0),
          sigma2    = stats::var(y1) / length(y1) + stats::var(y0) / length(y0)
        )
      }
    }
    est <- do.call(rbind, rows)

    # Per-stratum NPMLE via GLmix (Y_i | mu_i ~ N(mu_i, sigma_i^2), mu_i ~ G_g)
    npmle_g <- vector("list", G)
    for (g in seq_len(G)) {
      eg <- est[est$stratum == g, ]
      npmle_g[[g]] <- REBayes::GLmix(x = eg$theta_hat, sigma = sqrt(eg$sigma2))
    }

    strata_obs <- est$stratum
    s2_obs     <- est$sigma2
    counts_g   <- as.numeric(table(factor(strata_obs, levels = seq_len(G))))

    return(list(
      data  = est,
      X_obs = est$theta_hat,
      G_hat = npmle_g,

      prior_sampler = function(n) {
        g   <- sample(seq_len(G), size = n, replace = TRUE, prob = counts_g)
        out <- numeric(n)
        for (gg in seq_len(G)) {
          sel <- which(g == gg)
          if (!length(sel)) next
          f   <- npmle_g[[gg]]
          out[sel] <- sample(f$x, size = length(sel), replace = TRUE, prob = f$y)
        }
        out
      },

      likelihood_sampler = function(theta) {
        n  <- length(theta)
        s2 <- sample(s2_obs, size = n, replace = TRUE)
        stats::rnorm(n, mean = theta, sd = sqrt(s2))
      },

      feature_transform = function(x) x
    ))
  }


  # ------------------------------------------------------------------
  # Wang et al. 2024 (neural-g) - Arizona Medicare CABG length-of-stay
  #   Source: Hilbe (2011), data("azdrg112", package = "COUNT"); n = 1798.
  #   Model:  Y_i | lambda_i ~ Poisson(lambda_i),  lambda_i ~ G   (eq. 21).
  #   EB:     neural-g (the paper's headline estimator) via the neuralG
  #           package, a PyTorch model wrapped through reticulate.
  #   GOTCHA: neuralG's Poisson branch is spelled "Possion" (sic); passing
  #           "Poisson" silently leaves the likelihood unset. Its Poisson
  #           support grid is pinned to [0, max(y)] with n_grid points.
  # ------------------------------------------------------------------
  if (key == "wang_arizona") {

    if (!requireNamespace("COUNT", quietly = TRUE))
      stop("Please install the 'COUNT' package.")
    if (!requireNamespace("neuralG", quietly = TRUE))
      stop("Install neuralG: remotes::install_github('shijiew97/neuralG').")

    py_require(c("torch", "torchvision"))

    data("azdrg112", package = "COUNT", envir = environment())
    if (!"los" %in% names(azdrg112))
      stop("Expected a 'los' (length-of-stay) column in azdrg112.")
    y_obs <- as.numeric(azdrg112$los)

    # Pin the torch RNG for a more reproducible fit (CPU; see notes below).
    try(reticulate::py_run_string("import torch; torch.manual_seed(1)"),
        silent = TRUE)

    # Poisson has no nuisance parameter -> param = 0. dist MUST be "Possion".
    fit <- neuralG::neural_g(Y = y_obs, n = length(y_obs),
                             dist = "Possion",
                             param = 0,
                             n_grid = 100,
                             num_it = 8000,
                             lr = 0.0003,
                             lrdecay = 1,
                             lr_power = 0.2,
                             c = 0.6,
                             verb = 0)
    w   <- fit$prob / sum(fit$prob)

    if (abs(sum(fit$support * w) - mean(y_obs)) > 0.5 * mean(y_obs))
      warning("Fitted Poisson prior mean is far from mean(y); check the fit.")

    G_hat <- list(support = fit$support, weights = w)

    return(list(
      data  = azdrg112,
      X_obs = y_obs,
      G_hat = G_hat,
      prior_sampler = function(n) {
        sample(G_hat$support, size = n, replace = TRUE, prob = G_hat$weights)
      },
      likelihood_sampler = function(theta) {
        stats::rpois(length(theta), lambda = theta)
      },
      feature_transform = function(x) x
    ))
  }

  # ------------------------------------------------------------------
  # Wang et al. 2024 (neural-g) - Framingham Heart Study SBP measurement error
  #   Source: Carroll et al. (2006), shipped as data("framingham"); n = 1615.
  #   Model:  y_ij = mu_i + eta_i, eta ~ N(0, sigma^2) (homogeneous, eq 24);
  #           average two SBP2 replicates -> y_i ~ N(mu_i, sigma^2/2); mu_i ~ G.
  #           sigma^2 from the within-person difference (eq 25).
  #   EB:     neural-g via the neuralG package (PyTorch through reticulate).
  #           dist = "Gaussian-framing" is the package's purpose-built FHS
  #           branch: Normal(mu, param) likelihood with the support grid pinned
  #           to [50, 263] (the SBP range). Use "Gaussian" for a data-driven
  #           [min(y), max(y)] grid instead.
  # ------------------------------------------------------------------
  if (key == "wang_framingham") {

    if (!requireNamespace("neuralG", quietly = TRUE))
      stop("Install neuralG: remotes::install_github('shijiew97/neuralG').")
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("Install reticulate package').")

    py_require(c("torch", "torchvision"))

    data("framingham", package = "ebcalibration", envir = environment())

    y1 <- framingham$SBP21
    y2 <- framingham$SBP22
    sigma2_eta <- 0.5 * stats::var(y1 - y2)   # single-measurement error var (eq 25)
    y_avg      <- 0.5 * (y1 + y2)             # working obs ~ N(mu, sigma^2/2)
    sigma2_obs <- sigma2_eta / 2              # variance of the average (eq 26)
    sd_obs     <- sqrt(sigma2_obs)

    # --- Standardize so neural-g sees a ~unit scale (the regime it's tuned for).
    mu_y <- mean(y_avg); sd_y <- stats::sd(y_avg)
    z    <- (y_avg - mu_y) / sd_y
    sd_z <- sd_obs / sd_y                      # noise SD on the z-scale

    try(reticulate::py_run_string("import torch; torch.manual_seed(1)"),
        silent = TRUE)

    fit <- neuralG::neural_g(Y = z, n = length(z),
                             dist = "Gaussian", param = sd_z,
                             n_grid = 100, num_it = 8000,
                             lr = 0.0003, lrdecay = 1, lr_power = 0.2,
                             c = 0.6, verb = 1)        # verb=1: watch the loss move

    w <- fit$prob; w[!is.finite(w)] <- 0
    if (sum(w) <= 0)
      stop("neural-g returned all-zero weights - the fit failed. Use REBayes::GLmix.")
    w <- w / sum(w)

    support <- fit$support * sd_y + mu_y       # map the grid back to the SBP scale
    G_hat   <- list(support = support, weights = w)

    if (abs(sum(support * w) - mu_y) > 0.1 * abs(mu_y))
      warning("Fitted Gaussian prior mean is far from mean(y); check the fit.")

    return(list(
      data  = framingham,
      X_obs = y_avg,
      G_hat = G_hat,
      prior_sampler = function(n) {
        sample(G_hat$support, size = n, replace = TRUE, prob = G_hat$weights)
      },
      likelihood_sampler = function(theta) {
        stats::rnorm(length(theta), mean = theta, sd = sd_obs)
      },
      feature_transform = function(x) x
    ))
  }


  # ------------------------------------------------------------------
  # 25a. Dicker & Zhao 2016 - Golub leukemia, ALL class (NPMLE, REBayes)
  # ------------------------------------------------------------------
  .golub_class_means <- function(group = c("ALL", "AML")) {
    group <- match.arg(group)
    if (!requireNamespace("golubEsets", quietly = TRUE))
      stop("Install 'golubEsets' from Bioconductor: BiocManager::install('golubEsets')")
    if (!requireNamespace("Biobase", quietly = TRUE))
      stop("Install the 'Biobase' package (Bioconductor).")

    data("Golub_Train", package = "golubEsets", envir = environment())
    es  <- Golub_Train                                  # 7129 genes x 38 subjects
    Xg  <- Biobase::exprs(es)
    lab <- as.character(Biobase::pData(es)$ALL.AML)     # "ALL" / "AML"

    s <- apply(Xg, 1L, stats::sd)                       # standardize each gene to var 1
    s[s == 0] <- 1
    Xs <- Xg / s

    in_k <- lab == group
    list(xbar = as.numeric(rowMeans(Xs[, in_k, drop = FALSE])),  # length 7129
         n_k  = sum(in_k))
  }

  if (key == "dickerzhao_golub_all") {

    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install the 'REBayes' package.")

    gm    <- .golub_class_means("ALL")
    x_obs <- gm$xbar
    sd_k  <- 1 / sqrt(gm$n_k)                          # N(mu, 1/n_k)

    fit   <- REBayes::GLmix(x = x_obs, sigma = sd_k)
    G_hat <- list(support = fit$x, weights = fit$y)

    return(list(
      data  = data.frame(xbar = x_obs),
      X_obs = x_obs,
      G_hat = G_hat,
      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },
      likelihood_sampler = function(theta) {
        stats::rnorm(length(theta), mean = theta, sd = sd_k)
      },
      feature_transform = function(x) x                 # Gaussian, known var: identity
    ))
  }

  # ------------------------------------------------------------------
  # 25b. Dicker & Zhao 2016 - Golub leukemia, AML class (NPMLE, REBayes)
  # ------------------------------------------------------------------
  if (key == "dickerzhao_golub_aml") {

    if (!requireNamespace("REBayes", quietly = TRUE))
      stop("Please install the 'REBayes' package.")

    gm    <- .golub_class_means("AML")
    x_obs <- gm$xbar
    sd_k  <- 1 / sqrt(gm$n_k)

    fit   <- REBayes::GLmix(x = x_obs, sigma = sd_k)
    G_hat <- list(support = fit$x, weights = fit$y)

    return(list(
      data  = data.frame(xbar = x_obs),
      X_obs = x_obs,
      G_hat = G_hat,
      prior_sampler = function(n) {
        sample(G_hat$support, size = n,
               replace = TRUE, prob = G_hat$weights / sum(G_hat$weights))
      },
      likelihood_sampler = function(theta) {
        stats::rnorm(length(theta), mean = theta, sd = sd_k)
      },
      feature_transform = function(x) x
    ))
  }


  # ------------------------------------------------------------------
  # Fallback
  # ------------------------------------------------------------------
  stop(
    sprintf(
      "Unknown paper/dataset combination: '%s' / '%s'.\n",
      paper, dataset
    ),
    "Supported combinations:\n",
    "  get_data('Willwerscheid', 'wOBA_normal' | 'wOBA_unimodal' | 'wOBA_npmle')\n",
    "  get_data('Koenker', 'Norberg')\n",
    "  get_data('Koenker', 'Flies')\n",
    "  get_data('Narasimhan', 'Surgery')\n",
    "  get_data('Narasimhan', 'Shakespeare')\n",
    "  get_data('Gu', 'baseball')\n",
    "  get_data('Baek', 'brainvessel' | 'hardcandy' | 'deathnotice' |\n",
    "                   'traffic' | 'microbial' | 'shakespeare')\n",
    "  get_data('Jana', 'hockey_npmle' | 'hockey_hellinger' | 'hockey_chisquare')\n",
    "  get_data('You', 'STAR_gaussian' | 'STAR_npmle')\n",
    "  get_data('Wang', 'arizona')\n",
    "  get_data('Wang', 'framingham')\n",
    "  get_data('DickerZhao', 'golub_all')\n",
    "  get_data('DickerZhao', 'golub_aml')\n",
    call. = FALSE
  )
}
