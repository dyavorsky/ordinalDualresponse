# Joint likelihood of the ordinal dual-response model and the aggregate MLE.

par_to_cut <- function(par, P, W) {
  if (W == 2) return(par[P + 1])
  par[P + 1] + c(0, cumsum(exp(par[(P + 2):(P + W - 1)])))
}

cut_to_par <- function(cut) {
  if (length(cut) == 1) return(cut)
  c(cut[1], log(diff(cut)))
}

# Tail-stable interval probability G(hi) - G(lo) from standardized bounds.
ord_prob_z <- function(lo, hi, model) {
  if (model == "B") {
    stats::plogis(hi) * stats::plogis(-lo) * (-expm1(lo - hi))
  } else {
    ea <- exp(-lo)
    eb <- exp(-hi)
    p <- exp(-eb) * (-expm1(-(ea - eb)))
    # deep left tail: exp(-hi) overflows and both CDFs underflow to zero
    p[!is.finite(eb)] <- 0
    p
  }
}

ord_prob <- function(mubar, cut, y, model) {
  caug <- c(-Inf, cut, Inf)
  ord_prob_z(caug[y] - mubar, caug[y + 1L] - mubar, model)
}

#' Joint log-likelihood of the ordinal dual-response model
#'
#' Per task: MNL(inside goods) x \[G(c_y - mubar) - G(c_\{y-1\} - mubar)\],
#' with mubar the log-sum inclusive value and G the standard Gumbel CDF
#' (model "A") or logistic CDF (model "B").
#'
#' @param beta coefficient vector
#' @param cut increasing utility-scale cut points (length W - 1)
#' @param dat a dataset from [dr_simulate()] (aggregate: pooled tasks)
#' @return total log-likelihood (scalar)
#' @export
dr_loglik <- function(beta, cut, dat) {
  -negloglik(c(beta, cut_to_par(cut)), dat)
}

negloglik <- function(par, dat) {
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  beta <- par[1:P]
  cut <- par_to_cut(par, P, W)
  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))
  ll_choice <- V[cbind(seq_len(n), dat$jstar)] - logS
  p_ord <- ord_prob(logS, cut, dat$y, dat$model)
  -(sum(ll_choice) + sum(log(pmax(p_ord, 1e-312))))
}

num_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(k) {
    h <- eps * max(1, abs(x[k]))
    xp <- x; xp[k] <- xp[k] + h
    xm <- x; xm[k] <- xm[k] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(1))
}

#' Aggregate maximum-likelihood estimation
#'
#' @param dat a dataset from [dr_simulate()]
#' @param start optional start values on the working scale
#' @return an object of class `dr_mle` with elements beta, se_beta, cut,
#'   se_cut (delta method), alpha/se_alpha for model "A", nll, observed
#'   information eigenvalues, gradient norm
#' @export
dr_mle <- function(dat, start = NULL) {
  design <- dat$design
  P <- design$P
  W <- dat$W
  if (is.null(start)) {
    freq <- tabulate(dat$y, nbins = W)
    cumq <- pmin(pmax(cumsum(freq)[1:(W - 1)] / sum(freq), 1e-4), 1 - 1e-4)
    ginv <- if (dat$model == "B") stats::qlogis else function(q) -log(-log(q))
    cut0 <- log(design$J) + ginv(cumq)
    if (W > 2) for (k in 2:(W - 1)) cut0[k] <- max(cut0[k], cut0[k - 1] + 1e-3)
    start <- c(rep(0, P), cut_to_par(cut0))
  }
  fn <- function(p) negloglik(p, dat)
  opt <- stats::optim(start, fn, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-12))
  H <- stats::optimHess(opt$par, fn)
  g <- num_grad(fn, opt$par)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  cut_hat <- par_to_cut(opt$par, P, W)
  Jc <- matrix(0, W - 1, length(opt$par))
  Jc[, P + 1] <- 1
  if (W > 2) {
    d <- exp(opt$par[(P + 2):(P + W - 1)])
    for (w in 2:(W - 1)) Jc[w, (P + 2):(P + w)] <- d[1:(w - 1)]
  }
  Vpar <- tryCatch(solve(H), error = function(e)
    matrix(NA, length(opt$par), length(opt$par)))
  out <- list(
    beta = stats::setNames(opt$par[1:P], colnames(design$X)),
    se_beta = sqrt(diag(Vpar)[1:P]),
    cut = cut_hat, se_cut = sqrt(diag(Jc %*% Vpar %*% t(Jc))),
    nll = opt$value, convergence = opt$convergence,
    grad_max = max(abs(g)), eigen = ev,
    par = opt$par, vcov_par = Vpar, model = dat$model, W = W, P = P
  )
  if (dat$model == "A") {
    alpha_hat <- exp(-exp(-cut_hat))
    out$alpha <- alpha_hat
    out$se_alpha <- out$se_cut * alpha_hat * (-log(alpha_hat))
  }
  class(out) <- "dr_mle"
  out
}

#' @export
print.dr_mle <- function(x, ...) {
  cat(sprintf("ordinal dual-response MLE (model %s, W = %d)\n", x$model, x$W))
  tab <- data.frame(estimate = c(x$beta, x$cut),
                    se = c(x$se_beta, x$se_cut),
                    row.names = c(names(x$beta), paste0("c", seq_along(x$cut))))
  print(round(tab, 4))
  cat(sprintf("logLik %.2f | convergence %d | max|grad| %.1e | min eig(H) %.1f\n",
              -x$nll, x$convergence, x$grad_max, min(x$eigen)))
  invisible(x)
}

#' @export
coef.dr_mle <- function(object, ...) c(object$beta, stats::setNames(object$cut, paste0("c", seq_along(object$cut))))

#' Predicted exceedance probabilities P(y >= k) for new tasks
#' @param object a `dr_mle` fit
#' @param design a design list (same attribute layout)
#' @param k scale threshold
#' @param ... unused
#' @export
predict.dr_mle <- function(object, design, k, ...) {
  V <- matrix(design$X %*% object$beta, nrow = design$n_tasks,
              ncol = design$J, byrow = TRUE)
  m <- row_max(V)
  mubar <- m + log(rowSums(exp(V - m)))
  G <- if (object$model == "B") stats::plogis else function(z) exp(-exp(-z))
  1 - G(object$cut[k - 1] - mubar)
}
