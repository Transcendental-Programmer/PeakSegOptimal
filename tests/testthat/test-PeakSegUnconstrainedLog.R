library(testthat)
context("PeakSegUnconstrainedLog")
library(PeakSegOptimal)

test_that("basic functionality works", {
  # Simple test case
  data.vec <- as.integer(c(1, 10, 14, 13))
  fit <- PeakSegUnconstrainedLog(data.vec, rep(1L, 4), 0)
  
  # With no penalty, we should get segments for each data point
  expect_equal(length(unique(fit$mean.vec[fit$mean.vec > 0])), 4)
})

test_that("huge penalty recovers no changes", {
  data.vec <- as.integer(c(1, 10, 14, 5))
  fit <- PeakSegUnconstrainedLog(data.vec, rep(1L, 4), 1e6)
  
  # With huge penalty, should get a single mean
  valid_means <- fit$mean.vec[fit$mean.vec > 0]
  expect_equal(length(unique(valid_means)), 1)
  exp.mean <- mean(data.vec)
  expect_equal(valid_means[1], exp.mean, tolerance = 1e-5)
})

test_that("unconstrained model can have up-up or down-down changes", {
  # Data that would benefit from consecutive up changes
  data.vec <- as.integer(c(1, 5, 15, 20))
  fit <- PeakSegUnconstrainedLog(data.vec, rep(1L, 4), 1)
  
  # Extract segments
  is_change <- fit$ends.vec > 0
  break_vec <- fit$ends.vec[is_change]
  if(length(break_vec) > 0) {
    first <- c(1, break_vec+1)
    last <- c(break_vec, length(data.vec))
    mean_vec <- fit$mean.vec[1:length(first)]
    
    # Check if there are consecutive increases or decreases
    has_consecutive_changes <- FALSE
    if(length(mean_vec) > 2) {
      changes <- diff(mean_vec)
      has_consecutive_changes <- any(changes[1:(length(changes)-1)] * changes[2:length(changes)] > 0)
      expect_true(has_consecutive_changes)
    }
  }
})

test_that("unconstrained model vs constrained model", {
  # Create data that would be modeled differently with constraints
  data.vec <- as.integer(c(1, 10, 15, 10, 15, 5))
  
  # Compare unconstrained and constrained models
  unconstrained <- PeakSegUnconstrainedLog(data.vec, rep(1L, 6), 5)
  constrained <- PeakSegFPOP(data.vec, rep(1L, 6), 5)
  
  # Extract info from both models
  unc_cost <- min(unconstrained$cost.mat)
  con_cost <- min(constrained$cost.mat)
  
  # The unconstrained model should have lower or equal cost
  expect_true(unc_cost <= con_cost + 1e-6)
})

test_that("edge cases work properly", {
  # All zeros
  data.vec <- as.integer(rep(0, 5))
  expect_error({
    PeakSegUnconstrainedLog(data.vec, rep(1L, 5), 0)
  }, "data\\[i\\]=0 for all i")
  
  # Single value (all the same)
  data.vec <- as.integer(rep(5, 5))
  fit <- PeakSegUnconstrainedLog(data.vec, rep(1L, 5), 0)
  valid_means <- fit$mean.vec[fit$mean.vec > 0]
  expect_equal(length(unique(valid_means)), 1)
  expect_equal(valid_means[1], 5)
  
  # One different value
  data.vec <- as.integer(c(rep(5, 4), 10))
  fit <- PeakSegUnconstrainedLog(data.vec, rep(1L, 5), 0.1)
  valid_means <- fit$mean.vec[fit$mean.vec > 0]
  expect_true(length(unique(valid_means)) > 1)
})

test_that("weighted loss works as expected", {
  data.vec <- as.integer(c(1, 10, 14, 5))
  
  # With equal weights
  fit1 <- PeakSegUnconstrainedLog(data.vec, rep(1L, 4), 5)
  
  # With custom weights
  weight.vec <- c(1, 2, 1, 1)
  fit2 <- PeakSegUnconstrainedLog(data.vec, weight.vec, 5)
  
  # The models should be different
  expect_false(identical(fit1$mean.vec, fit2$mean.vec))
})
