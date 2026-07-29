## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(ordinalDR)
N <- 80; T_tasks <- 10
beta_bar <- c(a2 = 0.8, a3 = -0.5, b2 = 0.4, b3 = 1.0, price = -0.9)
Sigma <- diag(c(0.6, 0.6, 0.6, 0.6, 0.3))
cut <- c(-1.0, 0.2, 1.2, 2.2)

des <- dr_panel_design(N, T_tasks, J = 4, seed = 1)
B <- dr_draw_betas(N, beta_bar, Sigma, seed = 2)
dat <- dr_simulate(des, B, cut, model = "B", seed = 3)

fit <- dr_hb(dat, mcmc = list(R = 4000, burn = 2000, thin = 4),
             seed = 4, verbose = FALSE)
summary(fit, truth = c(beta_bar, cut))

## -----------------------------------------------------------------------------
bi <- dr_hb_betai(fit)
cor(as.numeric(bi), as.numeric(B))

sp <- dr_split_holdout(dat, 2)
fit_est <- dr_hb(sp$est, mcmc = list(R = 3000, burn = 1500, thin = 3),
                 seed = 5, verbose = FALSE)
dr_holdout_ll(fit_est, sp$holdout)$total
dr_waic(fit_est$ll_resp)$waic

## -----------------------------------------------------------------------------
fit_het <- dr_hb(dat, mcmc = list(R = 4000, burn = 2000, thin = 4),
                 het_cut = TRUE, seed = 6, verbose = FALSE)
round(colMeans(fit_het$draws$phibar), 3)

