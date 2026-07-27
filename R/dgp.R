# Simulation of the ordinal dual-response behavioral process.
# Data are always simulated from raw utilities (Gumbel draws, argmax, and the
# reporting rules), never from the derived likelihood, so that estimation
# validates the derivation as well as the code.

rgumbel <- function(n) -log(-log(runif(n)))

row_max <- function(M) do.call(pmax, as.data.frame(M))

#' Random attribute design for dual-response tasks
#'
#' Two 3-level attributes (dummy-coded against omitted reference levels) plus
#' a continuous price-like attribute; no intercept and no full dummy set, so
#' the identification condition rank(\[X 1\]) = P + 1 holds by construction.
#'
#' @param n_tasks number of tasks
#' @param J alternatives (inside goods) per task
#' @param seed RNG seed
#' @param intercept append a column of ones (deliberately violates
#'   identification; useful for demonstrations)
#' @return list with X, n_tasks, J, P
#' @export
dr_design <- function(n_tasks, J, seed, intercept = FALSE) {
  set.seed(seed)
  n_rows <- n_tasks * J
  a <- sample(1:3, n_rows, replace = TRUE)
  b <- sample(1:3, n_rows, replace = TRUE)
  price <- stats::runif(n_rows, 0.5, 2.5)
  X <- cbind(a2 = as.numeric(a == 2), a3 = as.numeric(a == 3),
             b2 = as.numeric(b == 2), b3 = as.numeric(b == 3), price = price)
  if (intercept) X <- cbind(X, const = 1)
  list(X = X, n_tasks = n_tasks, J = J, P = ncol(X))
}

#' Panel design: N respondents x T tasks
#'
#' @param N respondents
#' @param T_tasks tasks per respondent
#' @param J alternatives per task
#' @param seed RNG seed
#' @param intercept see [dr_design()]
#' @return a design list with respondent index vectors `resp` (per task) and
#'   `resp_row` (per design row)
#' @export
dr_panel_design <- function(N, T_tasks, J, seed, intercept = FALSE) {
  des <- dr_design(N * T_tasks, J, seed = seed, intercept = intercept)
  des$N <- N
  des$T_tasks <- T_tasks
  des$resp <- rep(seq_len(N), each = T_tasks)
  des$resp_row <- rep(des$resp, each = J)
  des
}

#' Check the identification rank condition rank(\[X 1\]) = P + 1
#' @param X stacked inside-good design matrix
#' @return list(rank, required, pass)
#' @export
dr_rank_check <- function(X) {
  aug <- cbind(X, 1)
  rk <- qr(aug)$rank
  list(rank = rk, required = ncol(aug), pass = rk == ncol(aug))
}

#' Draw individual coefficients from a multivariate-normal population
#' @param N respondents
#' @param beta_bar population mean
#' @param Sigma population covariance
#' @param seed RNG seed
#' @export
dr_draw_betas <- function(N, beta_bar, Sigma, seed) {
  set.seed(seed)
  P <- length(beta_bar)
  Z <- matrix(stats::rnorm(N * P), N, P)
  sweep(Z %*% chol(Sigma), 2, beta_bar, "+")
}

#' Simulate the ordinal dual-response behavioral process
#'
#' @param design a design from [dr_panel_design()] (use N = 1 for pooled
#'   aggregate data)
#' @param Bmat N x P matrix of individual coefficients
#' @param cut common cut-point vector (length W - 1, utility scale) or an
#'   N x (W - 1) matrix for heterogeneous cut points. For Model A, utility
#'   cuts map to probability cuts by alpha = exp(-exp(-cut)).
#' @param model "A" (probability report, Gumbel link) or "B" (graded
#'   comparison, logistic link)
#' @param seed RNG seed
#' @param model2 "inclusive" (the paper's second stage) or "2max" (DR-2Max
#'   style: latent is the chosen option's deterministic utility plus fresh
#'   logistic noise)
#' @param lambda slope on the inclusive value in the second-stage latent
#'   (Model B semantics; 1 is the theoretical value)
#' @param shock "gumbel" or "normal" (variance-matched; for misspecification
#'   studies)
#' @param drift NULL, `list(type = "linear", delta =)` or
#'   `list(type = "update", rho =)` for task-order drift in the effective cut
#'   points
#' @return a data list consumable by [dr_mle()] and [dr_hb()]
#' @export
dr_simulate <- function(design, Bmat, cut, model = c("A", "B"), seed,
                        model2 = c("inclusive", "2max"), lambda = 1,
                        shock = c("gumbel", "normal"), drift = NULL) {
  model <- match.arg(model)
  model2 <- match.arg(model2)
  shock <- match.arg(shock)
  set.seed(seed)
  n <- design$n_tasks
  J <- design$J
  N <- if (!is.null(design$N)) design$N else 1L
  if (is.null(design$resp)) {
    design$resp <- rep(1L, n)
    design$resp_row <- rep(1L, n * J)
    design$N <- 1L
    design$T_tasks <- n
  }
  stopifnot(nrow(Bmat) == N, ncol(Bmat) == design$P)

  cutmat <- if (is.matrix(cut)) cut else matrix(cut, N, length(cut), byrow = TRUE)
  W <- ncol(cutmat) + 1L

  Bx <- Bmat[design$resp_row, , drop = FALSE]
  V <- matrix(rowSums(design$X * Bx), nrow = n, ncol = J, byrow = TRUE)
  eta <- if (shock == "gumbel") rgumbel(n * J) else stats::rnorm(n * J, sd = pi / sqrt(6))
  u <- V + matrix(eta, n, J)
  jstar <- max.col(u, ties.method = "first")
  ustar <- u[cbind(seq_len(n), jstar)]

  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))

  latent <- if (model2 == "2max") {
    V[cbind(seq_len(n), jstar)] + stats::rlogis(n)
  } else if (lambda != 1) {
    stats::rlogis(n, location = lambda * logS)
  } else if (model == "A") {
    ustar
  } else {
    ustar - rgumbel(n)
  }

  shift <- numeric(n)
  if (!is.null(drift)) {
    T_tasks <- design$T_tasks
    tt <- rep(seq_len(T_tasks), times = design$N)
    if (drift$type == "linear") {
      shift <- drift$delta * (tt - 1)
    } else if (drift$type == "update") {
      psi <- numeric(n)
      for (i in seq_len(design$N)) {
        rows <- which(design$resp == i)
        for (k in seq_along(rows)[-1]) {
          psi[rows[k]] <- if (stats::runif(1) < stats::plogis(drift$rho))
            logS[rows[k - 1]] else psi[rows[k - 1]]
        }
      }
      shift <- psi
    }
  }

  eff_cut <- cutmat[design$resp, , drop = FALSE] + shift
  y <- 1L + rowSums(latent >= eff_cut)

  list(design = design, resp = design$resp, resp_row = design$resp_row,
       N = design$N, jstar = jstar, y = as.integer(y), W = W, model = model,
       Bmat_true = Bmat, cut_true = cut, model2 = model2, lambda = lambda,
       shock = shock, drift = drift)
}

#' Split a panel dataset into estimation and holdout tasks per respondent
#' @param dat a dataset from [dr_simulate()] (panel form)
#' @param n_holdout tasks per respondent to hold out (the last ones)
#' @export
dr_split_holdout <- function(dat, n_holdout) {
  des <- dat$design
  T_tasks <- des$T_tasks
  keep_t <- rep(seq_len(T_tasks) <= T_tasks - n_holdout, times = des$N)
  subset_panel <- function(keep) {
    rows <- rep(keep, each = des$J)
    d <- des
    d$X <- des$X[rows, , drop = FALSE]
    d$n_tasks <- sum(keep)
    d$T_tasks <- d$n_tasks / des$N
    d$resp <- des$resp[keep]
    d$resp_row <- rep(d$resp, each = des$J)
    out <- dat
    out$design <- d
    out$resp <- d$resp
    out$resp_row <- d$resp_row
    out$jstar <- dat$jstar[keep]
    out$y <- dat$y[keep]
    out
  }
  list(est = subset_panel(keep_t), holdout = subset_panel(!keep_t))
}

#' Dichotomize an ordinal dual response at a cut ("y >= k counts as buy")
#' @param dat a dataset from [dr_simulate()]
#' @param k the dichotomization category (e.g. W for top box, W - 1 for
#'   top-2-box on a 5-point scale use k = 4)
#' @export
dr_dichotomize <- function(dat, k) {
  dat$y <- 1L + (dat$y >= k)
  dat$W <- 2L
  dat
}
