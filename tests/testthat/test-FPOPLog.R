library(testthat)
library(PeakSegOptimal)

test_that("PeakSegFPOPLog returns one segment with huge penalty", {
  # Use data with variation to avoid the "all values same" error
  data_vec <- c(5, 6, 5, 7, 5)
  weight_vec <- rep(1, length(data_vec))
  data_count <- length(data_vec)
  penalty <- 1e6
  
  # Use the R wrapper function
  result <- PeakSegFPOPLog(data_vec, weight_vec, penalty)
  
  # With huge penalty, a single segment is favored.
  # Check that the segment ends at the last data point
  expect_equal(result$ends.vec[1], data_count)  # 1-indexed in R wrapper
  
  # And the mean should be finite
  expect_true(is.finite(result$mean.vec[1]))
})

test_that("PeakSegFPOPLog returns multiple segments with zero penalty", {
  data_vec <- c(5, 8, 5, 9, 7, 5, 10)
  weight_vec <- rep(1, length(data_vec))
  data_count <- length(data_vec)
  penalty <- 0
  
  # Use the R wrapper function
  result <- PeakSegFPOPLog(data_vec, weight_vec, penalty)
  
  # With zero penalty we expect segmentation to record at least one change point
  # First segment should not cover all data
  expect_true(result$ends.vec[1] < data_count)
})
