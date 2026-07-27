test_that("aggregate MLE recovers parameters (seeded smoke test)", {
  beta <- c(0.8, -0.5, 0.4, 1.0, -0.9)
  cut <- c(-1, 0.2, 1.2, 2.2)
  for (model in c("A", "B")) {
    des <- dr_design(6000, 4, seed = 11)
    dat <- dr_simulate(des, matrix(beta, 1, 5), cut, model = model, seed = 12)
    fit <- dr_mle(dat)
    expect_equal(fit$convergence, 0)
    expect_lt(max(abs(fit$beta - beta)), 0.12)
    expect_lt(max(abs(fit$cut - cut)), 0.15)
    expect_gt(min(fit$eigen), 0)
  }
})

test_that("dr_loglik matches negloglik", {
  des <- dr_design(500, 4, seed = 3)
  dat <- dr_simulate(des, matrix(c(0.5, -0.2, 0.3, 0.6, -0.7), 1, 5),
                     c(-1, 0, 1, 2), model = "B", seed = 4)
  ll <- dr_loglik(c(0.5, -0.2, 0.3, 0.6, -0.7), c(-1, 0, 1, 2), dat)
  expect_equal(ll, -ordinalDR:::negloglik(c(0.5, -0.2, 0.3, 0.6, -0.7,
                                            ordinalDR:::cut_to_par(c(-1, 0, 1, 2))), dat),
               tolerance = 1e-12)
})

test_that("Fisher information matches observed information at the truth", {
  beta <- c(0.8, -0.5, 0.4, 1.0, -0.9)
  cut <- c(-1, 0.2, 1.2, 2.2)
  des <- dr_design(20000, 4, seed = 21)
  dat <- dr_simulate(des, matrix(beta, 1, 5), cut, model = "B", seed = 22)
  I_exp <- dr_fisher(des, beta, cut, "B")
  H_obs <- stats::optimHess(c(beta, ordinalDR:::cut_to_par(cut)),
                            function(p) ordinalDR:::negloglik(p, dat))
  # compare in the beta block on the natural scale (cut blocks differ by the
  # working-parameter transform)
  expect_lt(max(abs(I_exp[1:5, 1:5] - H_obs[1:5, 1:5]) / abs(I_exp[1:5, 1:5])), 0.06)
})

test_that("omega is monotone under refinement and bounded by the latent info", {
  cut5 <- c(-1, 0.2, 1.2, 2.2)
  mubar <- seq(-2, 3, length.out = 11)
  for (model in c("A", "B")) {
    om5 <- dr_omega(mubar, cut5, model)
    om3 <- dr_omega(mubar, cut5[c(1, 3)], model)   # coarsening of the 5-pt partition
    om2 <- dr_omega(mubar, cut5[3], model)         # coarsening of the 3-pt partition
    lim <- if (model == "A") 1 else 1 / 3
    expect_true(all(om3 <= om5 + 1e-12))
    expect_true(all(om2 <= om3 + 1e-12))
    expect_true(all(om5 <= lim + 1e-12))
  }
})

test_that("specification tests behave on null and 2max data", {
  beta <- c(0.8, -0.5, 0.4, 1.0, -0.9)
  cut <- c(-1, 0.2, 1.2, 2.2)
  des <- dr_design(4000, 4, seed = 31)
  dat0 <- dr_simulate(des, matrix(beta, 1, 5), cut, model = "B", seed = 32)
  tb0 <- dr_spec_tests(dat0)
  expect_gt(tb0$p_lr[1], 0.001)               # lambda = 1 not rejected
  expect_gt(tb0$p_lr[2], 0.001)               # delta = 0 not rejected
  expect_lt(tb0$p_lr[3], 1e-6)                # DR-2Max rejected
  dat1 <- dr_simulate(des, matrix(beta, 1, 5), cut, model = "B", seed = 33,
                      model2 = "2max")
  tb1 <- dr_spec_tests(dat1)
  expect_lt(tb1$p_lr[2], 1e-6)                # delta = 0 rejected
  expect_gt(tb1$p_lr[3], 0.001)               # delta = 1 not rejected
  expect_lt(abs(tb1$estimate[2] - 1), 0.2)    # delta_hat near 1
})

test_that("HB sampler runs and produces sane output (tiny smoke)", {
  des <- dr_panel_design(30, 6, 4, seed = 41)
  B <- dr_draw_betas(30, c(0.8, -0.5, 0.4, 1.0, -0.9), diag(0.4, 5), seed = 42)
  dat <- dr_simulate(des, B, c(-1, 0.2, 1.2, 2.2), model = "B", seed = 43)
  fit <- dr_hb(dat, mcmc = list(R = 800, burn = 400, thin = 4),
               seed = 44, verbose = FALSE)
  expect_true(all(is.finite(fit$draws$phibar)))
  expect_true(all(is.finite(fit$draws$cut)))
  expect_gt(fit$accept$phi, 0.02)
  expect_lt(fit$accept$phi, 0.8)
  sp <- dr_split_holdout(dat, 1)
  fit2 <- dr_hb(sp$est, mcmc = list(R = 600, burn = 300, thin = 3),
                seed = 45, verbose = FALSE)
  h <- dr_holdout_ll(fit2, sp$holdout)
  expect_true(is.finite(h$total))
  w <- dr_waic(fit2$ll_resp)
  expect_true(is.finite(w$waic))
})

test_that("Nash solver finds the symmetric equilibrium of a symmetric market", {
  # simple symmetric logit market with outside margin: 3 identical goods
  shelf_fn <- function(p) cbind(price = p)
  BigB <- matrix(-1.5, 200, 1)
  dm <- function(p) dr_demand_dual(shelf_fn, p, BigB, cut_row = 0)
  eq <- dr_nash_prices(dm, cost = rep(0.5, 3))
  expect_lt(eq$spread, 0.02)
  expect_lt(max(abs(eq$p - mean(eq$p))), 0.01)           # symmetry
  expect_true(all(eq$p > 0.5))                           # positive markups
})
