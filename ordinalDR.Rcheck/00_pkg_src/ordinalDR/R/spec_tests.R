# Specification tests: the unit inclusive-value slope (lambda = 1) and the
# conditioning parameter spanning our model (delta = 0) and DR-2Max (delta = 1).

negloglik_ext <- function(par, dat, ext = c("lambda", "delta"), fix = NULL) {
  ext <- match.arg(ext)
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  if (is.null(fix)) {
    theta <- par[P + W]
    par_core <- par[1:(P + W - 1)]
  } else {
    theta <- fix
    par_core <- par
  }
  beta <- par_core[1:P]
  cut <- par_to_cut(par_core, P, W)
  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))
  Vstar <- V[cbind(seq_len(n), dat$jstar)]
  ll_choice <- Vstar - logS
  pred <- if (ext == "lambda") theta * logS else logS + theta * (Vstar - logS)
  p_ord <- ord_prob(pred, cut, dat$y, dat$model)
  -(sum(ll_choice) + sum(log(pmax(p_ord, 1e-312))))
}

#' Specification-test battery for the ordinal dual-response model
#'
#' Runs three one-degree-of-freedom tests from the aggregate likelihood:
#' the unit-slope test (H0: lambda = 1, an overidentification test of the
#' shared-utility account), the conditioning test (H0: delta = 0, the
#' inclusive-value second stage), and the reverse conditioning test
#' (H0: delta = 1, the DR-2Max second stage). Extended fits are warm-started
#' at the restricted MLE.
#'
#' @param dat aggregate dataset
#' @param mle optional precomputed [dr_mle()] fit
#' @return data frame with estimates, standard errors, Wald and LR statistics
#'   and p-values
#' @export
dr_spec_tests <- function(dat, mle = NULL) {
  if (is.null(mle)) mle <- dr_mle(dat)
  start_core <- mle$par
  fit_ext <- function(ext, theta0) {
    fn <- function(p) negloglik_ext(p, dat, ext = ext)
    opt <- stats::optim(c(start_core, theta0), fn, method = "BFGS",
                        control = list(maxit = 1000, reltol = 1e-12))
    H <- stats::optimHess(opt$par, fn)
    Vc <- tryCatch(solve(H), error = function(e)
      matrix(NA, length(opt$par), length(opt$par)))
    k <- length(opt$par)
    list(theta = opt$par[k], se = sqrt(Vc[k, k]), nll = opt$value)
  }
  fl <- fit_ext("lambda", 1)
  fd <- fit_ext("delta", 0)
  opt_d1 <- stats::optim(start_core,
                         function(p) negloglik_ext(p, dat, ext = "delta", fix = 1),
                         method = "BFGS",
                         control = list(maxit = 1000, reltol = 1e-12))
  lr <- c(2 * (mle$nll - fl$nll), 2 * (mle$nll - fd$nll),
          2 * (opt_d1$value - fd$nll))
  wz <- c((fl$theta - 1) / fl$se, fd$theta / fd$se, (fd$theta - 1) / fd$se)
  data.frame(
    test = c("lambda = 1", "delta = 0", "delta = 1 (DR-2Max)"),
    estimate = c(fl$theta, fd$theta, fd$theta),
    se = c(fl$se, fd$se, fd$se),
    wald_z = wz,
    lr = lr,
    p_wald = 2 * stats::pnorm(-abs(wz)),
    p_lr = stats::pchisq(lr, df = 1, lower.tail = FALSE)
  )
}
