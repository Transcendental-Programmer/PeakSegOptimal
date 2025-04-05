library(PeakSegOptimal)
library(microbenchmark)

# Create test data that would highlight differences in the approaches
set.seed(456)
n_points <- 100
x <- 1:n_points
peak_positions <- c(20, 50, 80)
peak_heights <- c(30, 50, 40)
baseline <- 5
signal <- baseline
# Create signal with multiple peaks
for(i in 1:length(peak_positions)) {
  signal <- signal + peak_heights[i] * exp(-(x - peak_positions[i])^2 / 50)
}
count.vec <- as.integer(rpois(n_points, signal))

# Run both algorithms
constrained <- PeakSegFPOPLog(count.vec, penalty=15)
unconstrained <- PeakSegUnconstrainedLog(count.vec, penalty=15)

# Compare results
cat("Cost comparison:\n")
cat("Constrained cost:", constrained$cost.mat[2*n_points], "\n")
cat("Unconstrained cost:", unconstrained$cost.mat[n_points], "\n")

# Benchmark performance
cat("\nPerformance comparison:\n")
benchmark <- microbenchmark(
  Constrained = PeakSegFPOPLog(count.vec, penalty=15),
  Unconstrained = PeakSegUnconstrainedLog(count.vec, penalty=15),
  times = 10
)
print(benchmark)
