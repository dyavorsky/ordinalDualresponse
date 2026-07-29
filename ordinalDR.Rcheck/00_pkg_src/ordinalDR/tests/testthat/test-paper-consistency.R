# Drift guard: the package must reproduce the paper repository's saved
# results. The paper's analysis code is standalone by design; this test
# regenerates one of its identification-test datasets with the package's
# generator (same seeds) and checks that the package MLE lands on the stored
# estimates.

test_that("package reproduces the paper repo's Model B set1 recovery fit", {
  paper <- Sys.getenv("ORDINALDR_PAPER_REPO",
                      file.path("..", "..", "..", "likert_dualresponse"))
  rds <- file.path(paper, "R", "output", "identification_test.rds")
  skip_if_not(file.exists(rds), "paper repo results not available")
  stored <- readRDS(rds)$recovery$B_set1_moderate
  skip_if(is.null(stored), "stored recovery entry not found")

  # settings mirrored from likert_dualresponse/R/identification_test.R:
  # Model B, set1_moderate is the 4th (model, set) combination -> seed_counter 4
  beta <- c(a2 = 0.8, a3 = -0.5, b2 = 0.4, b3 = 1.0, price = -0.9)
  cutB <- c(-1.0, 0.2, 1.2, 2.2)
  des <- dr_design(20000, 4, seed = 104)
  dat <- dr_simulate(des, matrix(beta, 1, 5), cutB, model = "B", seed = 504)
  fit <- dr_mle(dat)
  est_stored <- stored$table$est[1:9]
  expect_equal(unname(c(fit$beta, fit$cut)), est_stored, tolerance = 2e-3)
})
