# Model-comparison utilities: Newton-Raftery, WAIC, holdout likelihood,
# aggregate posterior sampling + bridge/Laplace marginal likelihoods.

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

#' Newton-Raftery log marginal density (harmonic mean), optionally trimmed
#' @param ll_draws total log-likelihood at each posterior draw
#' @param trim fraction of the smallest-likelihood draws to drop
#' @export
dr_lmd_nr <- function(ll_draws, trim = 0) {
  if (trim > 0) ll_draws <- ll_draws[ll_draws >= stats::quantile(ll_draws, trim)]
  -logmeanexp(-ll_draws)
}

#' WAIC from per-respondent pointwise log-likelihood draws
#' @param ll_resp_draws (ndraws x N) matrix
#' @export
dr_waic <- function(ll_resp_draws) {
  lppd <- sum(apply(ll_resp_draws, 2, logmeanexp))
  p_waic <- sum(apply(ll_resp_draws, 2, stats::var))
  list(waic = -2 * (lppd - p_waic), lppd = lppd, p_waic = p_waic,
       elpd = lppd - p_waic)
}

#' Log posterior-predictive density of holdout tasks
#' @param fit a `dr_hb` fit with stored individual draws
#' @param dat_holdout holdout dataset (same respondents; see
#'   [dr_split_holdout()])
#' @export
dr_holdout_ll <- function(fit, dat_holdout) {
  stopifnot(!is.null(fit$draws$Phi))
  P <- fit$settings$P
  W <- fit$settings$W
  het_cut <- fit$settings$het_cut
  nkeep <- dim(fit$draws$Phi)[1]
  N <- dim(fit$draws$Phi)[2]
  ll <- matrix(NA_real_, nkeep, N)
  for (r in seq_len(nkeep)) {
    Phi_r <- fit$draws$Phi[r, , , drop = TRUE]
    if (is.null(dim(Phi_r))) Phi_r <- matrix(Phi_r, N, 1)
    cutmat <- if (het_cut) phi_to_cutmat(Phi_r, P, W)
              else matrix(fit$draws$cut[r, ], N, W - 1, byrow = TRUE)
    ll[r, ] <- resp_loglik_all(dat_holdout, Phi_r[, 1:P, drop = FALSE], cutmat,
                               fit$settings$choice_only)
  }
  list(total = sum(apply(ll, 2, logmeanexp)),
       by_resp = apply(ll, 2, logmeanexp))
}

#' Posterior sampler for the aggregate model
#'
#' Random-walk Metropolis with the MLE's inverse Hessian as the proposal
#' shape and a N(0, tau2 I) prior on the working parameters. Used for exact
#' marginal-likelihood computation via bridge sampling.
#'
#' @param dat aggregate dataset
#' @param iter,burn,thin MCMC controls
#' @param tau2 prior variance
#' @param seed optional seed
#' @export
dr_bayes_agg <- function(dat, iter = 20000, burn = 5000, thin = 5,
                         tau2 = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mle <- dr_mle(dat)
  d <- length(mle$par)
  H <- stats::optimHess(mle$par, function(p) negloglik(p, dat))
  Vprop <- chol2inv(chol(H))
  Lp <- t(chol(Vprop)) * (2.4 / sqrt(d))
  lpost <- function(p) -negloglik(p, dat) - 0.5 * sum(p^2) / tau2
  cur <- mle$par
  cur_lp <- lpost(cur)
  nkeep <- floor((iter - burn) / thin)
  draws <- matrix(NA_real_, nkeep, d)
  acc <- 0; ki <- 0; sc <- 1
  for (it in seq_len(iter)) {
    prop <- cur + as.numeric(Lp %*% stats::rnorm(d)) * sc
    prop_lp <- lpost(prop)
    if (log(stats::runif(1)) < prop_lp - cur_lp) {
      cur <- prop; cur_lp <- prop_lp; acc <- acc + 1
    }
    if (it <= burn && it %% 200 == 0) {
      sc <- min(max(sc * exp(1.5 * (acc / 200 - 0.234)), 0.1), 10)
      acc <- 0
    }
    if (it == burn) acc <- 0
    if (it > burn && (it - burn) %% thin == 0) {
      ki <- ki + 1
      draws[ki, ] <- cur
    }
  }
  list(draws = draws, accept = acc / (iter - burn), mle = mle, tau2 = tau2,
       dat = dat)
}

#' Bridge-sampling log marginal likelihood for the aggregate model
#' @param fit_agg output of [dr_bayes_agg()]
#' @export
dr_bridge_lmd <- function(fit_agg) {
  if (!requireNamespace("bridgesampling", quietly = TRUE)) {
    warning("bridgesampling not installed; returning NA")
    return(NA_real_)
  }
  tau2 <- fit_agg$tau2
  lpost_fn <- function(pars, data) {
    -negloglik(pars, data) - 0.5 * sum(pars^2) / tau2 -
      0.5 * length(pars) * log(2 * pi * tau2)
  }
  samples <- fit_agg$draws
  colnames(samples) <- paste0("p", seq_len(ncol(samples)))
  lb <- rep(-Inf, ncol(samples)); ub <- rep(Inf, ncol(samples))
  names(lb) <- names(ub) <- colnames(samples)
  bs <- bridgesampling::bridge_sampler(samples = samples,
                                       log_posterior = lpost_fn,
                                       data = fit_agg$dat,
                                       lb = lb, ub = ub, silent = TRUE)
  bs$logml
}

#' Laplace-approximation log marginal likelihood for the aggregate model
#' @param fit_agg output of [dr_bayes_agg()]
#' @export
dr_laplace_lmd <- function(fit_agg) {
  mle <- fit_agg$mle
  d <- length(mle$par)
  H <- stats::optimHess(mle$par, function(p) negloglik(p, fit_agg$dat))
  lprior <- -0.5 * sum(mle$par^2) / fit_agg$tau2 -
    0.5 * d * log(2 * pi * fit_agg$tau2)
  as.numeric(-mle$nll + lprior + 0.5 * d * log(2 * pi) -
               0.5 * determinant(H, logarithm = TRUE)$modulus)
}
