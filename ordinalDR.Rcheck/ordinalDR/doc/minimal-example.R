## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(ordinalDR)
beta <- 1
cut <- c(0, 2)                       # W = 3: "unlikely" | "unsure" | "likely"
x1 <- c(2, 1)                        # attractive task
x2 <- c(-1, -2)                      # unattractive task

mnl <- function(x) exp(beta * x) / sum(exp(beta * x))
rbind(task1 = mnl(x1), task2 = mnl(x2))

## -----------------------------------------------------------------------------
mubar <- c(log(sum(exp(beta * x1))), log(sum(exp(beta * x2))))
probs <- function(mu, model) diff(c(0, if (model == "B") plogis(cut - mu)
                                    else exp(-exp(-(cut - mu))), 1))
rbind(task1_B = probs(mubar[1], "B"), task2_B = probs(mubar[2], "B"),
      task1_A = probs(mubar[1], "A"), task2_A = probs(mubar[2], "A"))

## -----------------------------------------------------------------------------
n <- 4000
X <- matrix(rep(c(x1, x2), n / 2), ncol = 1)
des <- list(X = X, n_tasks = n, J = 2, P = 1)
dat <- dr_simulate(des, matrix(beta, 1, 1), cut, model = "B", seed = 7)
fit <- dr_mle(dat)
fit

## -----------------------------------------------------------------------------
dr_rank_check(des$X)

