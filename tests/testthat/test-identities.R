test_that("interval probabilities sum to one for both links, including extreme cuts", {
  for (model in c("A", "B")) {
    for (cut in list(c(-1, 0.2, 1.2, 2.2), c(-8, -6, 7, 9), c(0.001, 0.002, 0.003, 0.004))) {
      mubar <- seq(-6, 6, length.out = 25)
      p <- vapply(1:5, function(w) ordinalDualresponse:::ord_prob(mubar, cut, rep(w, 25), model),
                  numeric(25))
      expect_true(all(p >= 0))
      expect_equal(rowSums(p), rep(1, 25), tolerance = 1e-12)
    }
  }
})

test_that("tail-stable interval probabilities have no NaNs at extremes", {
  for (model in c("A", "B")) {
    p <- ordinalDualresponse:::ord_prob_z(c(-Inf, -750, 700), c(-745, -740, Inf), model)
    expect_true(all(is.finite(p)))
    expect_true(all(p >= 0 & p <= 1))
  }
})

test_that("W = 2 Model B collapses exactly to MNL with a no-choice constant", {
  set.seed(1)
  V <- c(0.7, -0.3, 1.1)
  S <- sum(exp(V))
  gam <- 0.4
  # joint P(choose j, buy) from the model: MNL x (1 - Lambda(gamma - mubar))
  mubar <- log(S)
  p_buy_j <- exp(V) / S * (1 - plogis(gam - mubar))
  # telescoped standard MNL with none ASC gamma
  expect_equal(p_buy_j, exp(V) / (exp(gam) + S), tolerance = 1e-12)
  # and P(j, not buy) = MNL x Lambda(gamma - mubar)
  p_nobuy_j <- exp(V) / S * plogis(gam - mubar)
  expect_equal(p_nobuy_j, exp(V) / S * exp(gam) / (exp(gam) + S), tolerance = 1e-12)
})

test_that("cut-point transform round-trips", {
  cut <- c(-1.3, -0.2, 0.9, 2.5)
  expect_equal(ordinalDualresponse:::par_to_cut(c(rep(0, 3), ordinalDualresponse:::cut_to_par(cut)), 3, 5),
               cut, tolerance = 1e-12)
})

test_that("identification rank check flags intercepts", {
  des <- dr_design(200, 4, seed = 5)
  expect_true(dr_rank_check(des$X)$pass)
  des_bad <- dr_design(200, 4, seed = 5, intercept = TRUE)
  expect_false(dr_rank_check(des_bad$X)$pass)
})

test_that("simulated joint pmf matches the closed-form likelihood", {
  # one fixed task, many behavioral draws vs MNL x ordinal product
  beta <- c(0.8, -0.5, 0.4, 1.0, -0.9)
  cut <- c(-1, 0.2, 1.2, 2.2)
  X0 <- rbind(c(1, 0, 0, 0, 1.0), c(0, 1, 1, 0, 2.0),
              c(0, 0, 0, 1, 0.8), c(0, 0, 1, 0, 1.6))
  des <- list(X = do.call(rbind, replicate(20000, X0, simplify = FALSE)),
              n_tasks = 20000, J = 4, P = 5)
  for (model in c("A", "B")) {
    dat <- dr_simulate(des, matrix(beta, 1, 5), cut, model = model, seed = 42)
    V <- as.numeric(X0 %*% beta)
    mnl <- exp(V) / sum(exp(V))
    mubar <- log(sum(exp(V)))
    pw <- ordinalDualresponse:::ord_prob(rep(mubar, 5), cut, 1:5, model)
    emp <- as.numeric(table(factor(dat$jstar, 1:4), factor(dat$y, 1:5))) / 20000
    theo <- as.numeric(outer(mnl, pw))
    expect_lt(max(abs(emp - theo)), 0.012)
  }
})
