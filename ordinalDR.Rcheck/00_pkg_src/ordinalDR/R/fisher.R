# Exact Fisher information and scale-design tools.

link_cdf <- function(z, model) if (model == "B") stats::plogis(z) else exp(-exp(-z))
link_pdf <- function(z, model) {
  if (model == "B") stats::dlogis(z) else exp(-(z + exp(-z)))
}

fisher_task <- function(Xt, beta, cut, model) {
  P <- ncol(Xt)
  W <- length(cut) + 1L
  K <- P + W - 1L
  V <- as.numeric(Xt %*% beta)
  eV <- exp(V - max(V))
  p <- eV / sum(eV)
  mubar <- max(V) + log(sum(eV))
  I <- matrix(0, K, K)
  I[1:P, 1:P] <- t(Xt) %*% (diag(p) - tcrossprod(p)) %*% Xt
  z <- cut - mubar
  Gz <- c(0, link_cdf(z, model), 1)
  gz <- c(0, link_pdf(z, model), 0)
  pi_w <- diff(Gz)
  xbar <- as.numeric(t(p) %*% Xt)
  for (w in seq_len(W)) {
    s <- numeric(K)
    s[1:P] <- -(gz[w + 1] - gz[w]) / pi_w[w] * xbar
    if (w <= W - 1) s[P + w] <- s[P + w] + gz[w + 1] / pi_w[w]
    if (w >= 2)     s[P + w - 1] <- s[P + w - 1] - gz[w] / pi_w[w]
    I <- I + pi_w[w] * tcrossprod(s)
  }
  I
}

#' Total Fisher information of a dual-response design
#'
#' Sums the per-task information (MNL block plus the ordinal-stage block over
#' the full (beta, cut) parameter) over all tasks in the design.
#'
#' @param design a design list
#' @param beta,cut evaluation parameters
#' @param model link ("A" or "B")
#' @return (P + W - 1) square information matrix
#' @export
dr_fisher <- function(design, beta, cut, model) {
  n <- design$n_tasks
  J <- design$J
  K <- design$P + length(cut)
  I <- matrix(0, K, K)
  for (t in seq_len(n)) {
    rows <- ((t - 1) * J + 1):(t * J)
    I <- I + fisher_task(design$X[rows, , drop = FALSE], beta, cut, model)
  }
  I
}

#' Second-stage information weight omega(mubar; cut)
#'
#' The scalar in the paper's information decomposition: the ordinal stage
#' contributes omega * xbar xbar' of information about beta per task.
#'
#' @param mubar vector of inclusive values
#' @param cut cut points
#' @param model link
#' @export
dr_omega <- function(mubar, cut, model) {
  z <- outer(mubar, cut, function(m, c) c - m)
  Gz <- cbind(0, link_cdf(z, model), 1)
  gz <- cbind(0, link_pdf(z, model), 0)
  W <- length(cut) + 1L
  om <- 0
  for (w in seq_len(W)) {
    om <- om + (gz[, w + 1] - gz[, w])^2 / (Gz[, w + 1] - Gz[, w])
  }
  om
}

mean_omega <- function(cut, model, m0, s, ngrid = 81) {
  z <- seq(-4, 4, length.out = ngrid)
  w <- stats::dnorm(z); w <- w / sum(w)
  sum(w * dr_omega(m0 + s * z, cut, model))
}

#' Cut points with equal expected category shares
#'
#' Population quantiles of the second-stage latent when the inclusive value is
#' N(m0, s^2) across tasks.
#'
#' @param W scale points
#' @param model link
#' @param m0,s inclusive-value distribution
#' @export
dr_equal_share_cuts <- function(W, model, m0, s) {
  qs <- seq_len(W - 1) / W
  Fm <- function(c) stats::integrate(function(m) link_cdf(c - m, model) *
                                       stats::dnorm(m, m0, s),
                                     m0 - 6 * s - 10, m0 + 6 * s + 10)$value
  vapply(qs, function(q) stats::uniroot(function(c) Fm(c) - q,
                                        c(m0 - 20, m0 + 20))$root, numeric(1))
}

#' Information-optimal cut points
#'
#' Maximizes the mean second-stage information weight over the inclusive-value
#' distribution.
#'
#' @inheritParams dr_equal_share_cuts
#' @return list(cut, omega)
#' @export
dr_optimal_cuts <- function(W, model, m0, s) {
  start <- dr_equal_share_cuts(W, model, m0, s)
  if (W == 2) {
    opt <- stats::optimize(function(c) -mean_omega(c, model, m0, s),
                           interval = c(m0 - 8, m0 + 8))
    return(list(cut = opt$minimum, omega = -opt$objective))
  }
  par0 <- c(start[1], log(diff(start)))
  fn <- function(p) -mean_omega(p[1] + c(0, cumsum(exp(p[-1]))), model, m0, s)
  opt <- stats::optim(par0, fn, method = "BFGS", control = list(maxit = 500))
  list(cut = opt$par[1] + c(0, cumsum(exp(opt$par[-1]))), omega = -opt$value)
}

#' Sample size for a target standard error
#'
#' Scales the per-task Fisher information of a reference design to find the
#' number of tasks (or respondents, given tasks per respondent) delivering a
#' target asymptotic standard error for a chosen parameter.
#'
#' @param design reference design (its size sets the information rate)
#' @param beta,cut,model evaluation point
#' @param param index of the target parameter in (beta, cut)
#' @param target_se desired standard error
#' @param per_respondent tasks per respondent (returns respondents if given)
#' @return list with n_tasks (and n_respondents if applicable) and the
#'   reference-design standard error
#' @export
dr_power <- function(design, beta, cut, model, param, target_se,
                     per_respondent = NULL) {
  I <- dr_fisher(design, beta, cut, model)
  Vc <- solve(I)
  se_ref <- sqrt(Vc[param, param])
  n_tasks <- design$n_tasks * (se_ref / target_se)^2
  out <- list(se_reference = se_ref, n_tasks = ceiling(n_tasks))
  if (!is.null(per_respondent)) out$n_respondents <- ceiling(n_tasks / per_respondent)
  out
}
