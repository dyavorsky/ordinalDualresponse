# Nash-pricing counterfactual helpers.

#' Dual-response shelf demand from stacked coefficient draws
#'
#' Demand for each of the alternatives on a shelf: MNL share times the
#' probability the purchase evaluation clears the buy threshold, averaged over
#' the rows of `BigB` (stacked posterior draws x respondents, or true
#' individual coefficients).
#'
#' @param shelf_fn function(p) returning the shelf design matrix (J x P) at
#'   price vector p
#' @param p price vector
#' @param BigB stacked coefficient matrix
#' @param cut_row buy threshold per row (scalar recycled, or vector matching
#'   draws)
#' @param link "A" or "B"
#' @export
dr_demand_dual <- function(shelf_fn, p, BigB, cut_row, link = "B") {
  V <- BigB %*% t(shelf_fn(p))
  m <- row_max(V)
  eV <- exp(V - m)
  S <- rowSums(eV)
  mnl <- eV / S
  mubar <- m + log(S)
  G <- if (link == "B") stats::plogis else function(z) exp(-exp(-z))
  pbuy <- 1 - G(cut_row - mubar)
  colMeans(mnl * pbuy)
}

#' Constant-incidence heuristic demand (MNL shares x fixed incidence)
#' @inheritParams dr_demand_dual
#' @param incidence fixed purchase incidence applied to every task
#' @export
dr_demand_constant <- function(shelf_fn, p, BigB, incidence) {
  V <- BigB %*% t(shelf_fn(p))
  m <- row_max(V)
  eV <- exp(V - m)
  mnl <- eV / rowSums(eV)
  colMeans(mnl) * incidence
}

#' Nash equilibrium prices by damped best-response iteration
#'
#' @param demand_fn function(p) returning the demand vector
#' @param cost unit-cost vector
#' @param starts list of starting price vectors (multi-start check)
#' @param lower,upper price bounds for the per-firm best response
#' @param damp damping factor in (0, 1]
#' @param tol convergence tolerance on max price change
#' @param maxit maximum iterations
#' @return list(p, spread) -- equilibrium prices and the max disagreement
#'   across starts
#' @export
dr_nash_prices <- function(demand_fn, cost, starts = NULL,
                           lower = NULL, upper = NULL,
                           damp = 0.5, tol = 1e-4, maxit = 200) {
  nb <- length(cost)
  if (is.null(lower)) lower <- cost
  if (is.null(upper)) upper <- cost + 4
  if (is.null(starts)) starts <- list(cost + 0.3, cost + 1, cost + 2)
  eqs <- lapply(starts, function(p0) {
    p <- p0
    for (it in seq_len(maxit)) {
      p_new <- p
      for (b in seq_len(nb)) {
        br <- stats::optimize(function(pb) {
          pp <- p_new; pp[b] <- pb
          -(demand_fn(pp)[b] * (pb - cost[b]))
        }, interval = c(lower[b], upper[b]))
        p_new[b] <- (1 - damp) * p[b] + damp * br$minimum
      }
      if (max(abs(p_new - p)) < tol) break
      p <- p_new
    }
    p
  })
  spread <- max(vapply(eqs, function(e) max(abs(e - eqs[[1]])), numeric(1)))
  if (spread > 0.02) warning(sprintf("equilibrium multi-start spread %.3f", spread))
  list(p = eqs[[1]], spread = spread)
}
