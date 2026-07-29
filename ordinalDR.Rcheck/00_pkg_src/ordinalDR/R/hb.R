# Hierarchical Bayesian estimation: random-walk Metropolis within Gibbs.

phi_to_cutmat <- function(Phi, P, W) {
  c1 <- Phi[, P + 1]
  if (W == 2) return(matrix(c1, ncol = 1))
  D <- exp(Phi[, (P + 2):(P + W - 1), drop = FALSE])
  cutmat <- matrix(0, nrow(Phi), W - 1)
  cutmat[, 1] <- c1
  for (w in 2:(W - 1)) cutmat[, w] <- cutmat[, w - 1] + D[, w - 1]
  cutmat
}

resp_loglik_all <- function(dat, Bmat, cutmat, choice_only = FALSE) {
  des <- dat$design
  n <- des$n_tasks
  J <- des$J
  Bx <- Bmat[dat$resp_row, , drop = FALSE]
  V <- matrix(rowSums(des$X * Bx), nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))
  ll <- V[cbind(seq_len(n), dat$jstar)] - logS
  if (!choice_only) {
    cut_aug <- cbind(-Inf, cutmat, Inf)
    lo <- cut_aug[cbind(dat$resp, dat$y)] - logS
    hi <- cut_aug[cbind(dat$resp, dat$y + 1L)] - logS
    ll <- ll + log(pmax(ord_prob_z(lo, hi, dat$model), 1e-312))
  }
  as.numeric(rowsum(ll, dat$resp, reorder = TRUE))
}

#' Default hierarchical priors
#' @param dim_phi dimension of the individual-parameter vector
#' @export
dr_default_priors <- function(dim_phi) {
  list(phibar0 = rep(0, dim_phi), A = 1 / 100, nu = dim_phi + 3,
       V0 = (dim_phi + 3) * diag(dim_phi), zeta_prec = 1 / 100)
}

draw_phibar <- function(Phi, Sigma, priors) {
  N <- nrow(Phi)
  Sinv <- chol2inv(chol(Sigma))
  prec <- N * Sinv + priors$A * diag(ncol(Phi))
  Vb <- chol2inv(chol(prec))
  m <- Vb %*% (Sinv %*% (N * colMeans(Phi)) + priors$A * priors$phibar0)
  as.numeric(m + t(chol(Vb)) %*% stats::rnorm(ncol(Phi)))
}

draw_sigma <- function(Phi, phibar, priors) {
  R <- sweep(Phi, 2, phibar)
  Vn <- priors$V0 + crossprod(R)
  nun <- priors$nu + nrow(Phi)
  Winv <- chol2inv(chol(Vn))
  Wdraw <- stats::rWishart(1, df = nun, Sigma = Winv)[, , 1]
  chol2inv(chol(Wdraw))
}

#' Hierarchical Bayesian estimation of the ordinal dual-response model
#'
#' Random-walk Metropolis within Gibbs (Rossi, Allenby & McCulloch style):
#' all respondents' individual parameters updated simultaneously by
#' independent RW proposals with covariance s_i^2 Sigma; the common cut-point
#' block by RW-MH on the unconstrained scale; population moments by conjugate
#' normal / inverse-Wishart draws. Step sizes adapt during burn-in only.
#'
#' @param dat a panel dataset from [dr_simulate()]
#' @param mcmc list(R, burn, thin)
#' @param het_cut individual-specific cut points (stacked with beta_i in one
#'   MVN hierarchy)
#' @param choice_only drop the ordinal stage (hierarchical MNL; validation)
#' @param priors see [dr_default_priors()]
#' @param keep_phi store individual-parameter draws (needed for individual
#'   metrics and holdout prediction)
#' @param seed optional RNG seed
#' @param verbose progress messages
#' @return object of class `dr_hb`
#' @export
dr_hb <- function(dat, mcmc = list(R = 30000, burn = 10000, thin = 10),
                  het_cut = FALSE, choice_only = FALSE, priors = NULL,
                  keep_phi = TRUE, seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  des <- dat$design
  N <- dat$N
  P <- des$P
  W <- dat$W
  dim_phi <- if (het_cut) P + W - 1 else P
  if (is.null(priors)) priors <- dr_default_priors(dim_phi)

  R_iter <- mcmc$R
  burn <- mcmc$burn
  thin <- mcmc$thin
  nkeep <- floor((R_iter - burn) / thin)

  Phi <- matrix(0, N, dim_phi)
  freq <- tabulate(dat$y, nbins = W)
  cumq <- pmin(pmax(cumsum(freq)[1:(W - 1)] / sum(freq), 1e-4), 1 - 1e-4)
  ginv <- if (dat$model == "B") stats::qlogis else function(q) -log(-log(q))
  cut0 <- log(des$J) + ginv(cumq)
  if (W > 2) for (k in 2:(W - 1)) cut0[k] <- max(cut0[k], cut0[k - 1] + 1e-3)
  zeta <- cut_to_par(cut0)
  if (het_cut) Phi[, (P + 1):(P + W - 1)] <- matrix(zeta, N, W - 1, byrow = TRUE)

  phibar <- colMeans(Phi)
  Sigma <- diag(dim_phi)
  cur_cutmat <- if (het_cut) phi_to_cutmat(Phi, P, W)
                else matrix(par_to_cut(zeta, 0, W), N, W - 1, byrow = TRUE)
  cur_ll <- resp_loglik_all(dat, Phi[, 1:P, drop = FALSE], cur_cutmat, choice_only)

  s_phi <- rep(2.93 / sqrt(dim_phi), N)
  acc_phi <- rep(0, N)
  s_zeta <- 0.05
  acc_zeta <- 0
  window <- 100

  keep_betabar <- matrix(NA_real_, nkeep, dim_phi)
  keep_Sigma <- array(NA_real_, c(nkeep, dim_phi, dim_phi))
  keep_cut <- if (!het_cut && !choice_only) matrix(NA_real_, nkeep, W - 1) else NULL
  keep_ll <- numeric(nkeep)
  keep_ll_resp <- matrix(NA_real_, nkeep, N)
  keep_Phi <- if (keep_phi) array(NA_real_, c(nkeep, N, dim_phi)) else NULL

  t0 <- Sys.time()
  ki <- 0
  for (it in seq_len(R_iter)) {
    U <- chol(Sigma)
    Z <- matrix(stats::rnorm(N * dim_phi), N, dim_phi) %*% U
    Prop <- Phi + Z * s_phi
    prop_cutmat <- if (het_cut) phi_to_cutmat(Prop, P, W) else cur_cutmat
    prop_ll <- resp_loglik_all(dat, Prop[, 1:P, drop = FALSE], prop_cutmat, choice_only)
    Sinv <- chol2inv(chol(Sigma))
    dcur <- sweep(Phi, 2, phibar)
    dprop <- sweep(Prop, 2, phibar)
    qcur <- rowSums((dcur %*% Sinv) * dcur)
    qprop <- rowSums((dprop %*% Sinv) * dprop)
    accept <- log(stats::runif(N)) < (prop_ll - 0.5 * qprop) - (cur_ll - 0.5 * qcur)
    Phi[accept, ] <- Prop[accept, ]
    cur_ll[accept] <- prop_ll[accept]
    if (het_cut && any(accept)) cur_cutmat[accept, ] <- prop_cutmat[accept, , drop = FALSE]
    acc_phi <- acc_phi + accept

    if (!het_cut && !choice_only) {
      zeta_prop <- zeta + stats::rnorm(W - 1, sd = s_zeta)
      cut_prop <- matrix(par_to_cut(zeta_prop, 0, W), N, W - 1, byrow = TRUE)
      ll_prop <- resp_loglik_all(dat, Phi[, 1:P, drop = FALSE], cut_prop, choice_only)
      lp <- (sum(ll_prop) - 0.5 * priors$zeta_prec * sum(zeta_prop^2)) -
            (sum(cur_ll) - 0.5 * priors$zeta_prec * sum(zeta^2))
      if (log(stats::runif(1)) < lp) {
        zeta <- zeta_prop
        cur_cutmat <- cut_prop
        cur_ll <- ll_prop
        acc_zeta <- acc_zeta + 1
      }
    }

    phibar <- draw_phibar(Phi, Sigma, priors)
    Sigma <- draw_sigma(Phi, phibar, priors)

    if (it <= burn && it %% window == 0) {
      s_phi <- pmin(pmax(s_phi * exp(1.5 * (acc_phi / window - 0.23)), 0.01), 10)
      acc_phi <- rep(0, N)
      if (!het_cut && !choice_only) {
        s_zeta <- min(max(s_zeta * exp(1.5 * (acc_zeta / window - 0.30)), 1e-3), 2)
        acc_zeta <- 0
      }
    }
    if (it == burn) { acc_phi <- rep(0, N); acc_zeta <- 0 }

    if (it > burn && (it - burn) %% thin == 0) {
      ki <- ki + 1
      keep_betabar[ki, ] <- phibar
      keep_Sigma[ki, , ] <- Sigma
      if (!is.null(keep_cut)) keep_cut[ki, ] <- par_to_cut(zeta, 0, W)
      keep_ll[ki] <- sum(cur_ll)
      keep_ll_resp[ki, ] <- cur_ll
      if (keep_phi) keep_Phi[ki, , ] <- Phi
    }
    if (verbose && it %% 5000 == 0) {
      message(sprintf("  iter %d/%d (%.1f min)", it, R_iter,
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }

  post_iters <- R_iter - burn
  structure(list(
    draws = list(phibar = keep_betabar, Sigma = keep_Sigma, cut = keep_cut,
                 Phi = keep_Phi),
    ll = keep_ll, ll_resp = keep_ll_resp,
    accept = list(phi = mean(acc_phi / post_iters),
                  zeta = if (het_cut || choice_only) NA else acc_zeta / post_iters),
    settings = list(mcmc = mcmc, het_cut = het_cut, choice_only = choice_only,
                    P = P, W = W, N = N, model = dat$model, priors = priors),
    runtime_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ), class = "dr_hb")
}

#' Posterior-mean individual coefficients
#' @param fit a `dr_hb` fit (with keep_phi = TRUE)
#' @export
dr_hb_betai <- function(fit) {
  apply(fit$draws$Phi[, , 1:fit$settings$P, drop = FALSE], c(2, 3), mean)
}

#' @export
summary.dr_hb <- function(object, truth = NULL, ...) {
  pm <- colMeans(object$draws$phibar)
  psd <- apply(object$draws$phibar, 2, stats::sd)
  out <- data.frame(param = paste0("phibar", seq_along(pm)),
                    post_mean = pm, post_sd = psd)
  if (!is.null(object$draws$cut)) {
    out <- rbind(out, data.frame(
      param = paste0("c", seq_len(ncol(object$draws$cut))),
      post_mean = colMeans(object$draws$cut),
      post_sd = apply(object$draws$cut, 2, stats::sd)))
  }
  if (!is.null(truth)) {
    out$truth <- truth
    out$z <- (out$post_mean - truth) / out$post_sd
  }
  out
}

#' @export
print.dr_hb <- function(x, ...) {
  cat(sprintf("ordinal dual-response HB fit (model %s, W = %d, N = %d%s%s)\n",
              x$settings$model, x$settings$W, x$settings$N,
              if (x$settings$het_cut) ", heterogeneous cuts" else "",
              if (x$settings$choice_only) ", choice only" else ""))
  cat(sprintf("draws kept %d | accept: phi %.2f, zeta %s | runtime %.1f min\n",
              nrow(x$draws$phibar), x$accept$phi,
              ifelse(is.na(x$accept$zeta), "-", sprintf("%.2f", x$accept$zeta)),
              x$runtime_min))
  print(round(summary(x)[, c("param", "post_mean", "post_sd")], 4))
  invisible(x)
}
