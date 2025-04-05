library(PeakSegOptimal)
library(ggplot2)
library(data.table)

# Function to create synthetic data with peaks
generate_synthetic_data <- function(n_points, n_peaks, noise_level=0.2) {
  set.seed(123)
  x <- 1:n_points
  baseline <- 5
  signal <- rep(baseline, n_points)
  
  # Add peaks at random locations
  peak_positions <- round(seq(n_points*0.1, n_points*0.9, length.out=n_peaks))
  peak_heights <- runif(n_peaks, 20, 50)
  peak_widths <- runif(n_peaks, n_points/30, n_points/15)
  
  for(i in 1:n_peaks) {
    signal <- signal + peak_heights[i] * exp(-((x - peak_positions[i])^2) / peak_widths[i])
  }
  
  # Convert to count data
  count.vec <- rpois(n_points, signal)
  
  list(
    count.vec = count.vec,
    true_peaks = peak_positions
  )
}

# Compare unconstrained vs constrained models across penalty values
compare_models <- function(count.vec, penalty_range) {
  results <- data.table()
  
  for(pen in penalty_range) {
    # Run unconstrained model
    unconstrained <- PeakSegUnconstrainedLog(count.vec, penalty = pen)
    
    # Extract segment info
    segments_unc <- sum(unconstrained$mean.vec > 0, na.rm=TRUE)
    cost_unc <- unconstrained$cost.mat[length(count.vec)]
    
    # Store results
    results <- rbind(results, data.table(
      penalty = pen,
      model = "Unconstrained",
      segments = segments_unc,
      cost = cost_unc
    ))
    
    # Try to run constrained (FPOP) model if it's available in the namespace
    if(exists("PeakSegFPOPLog")) {
      constrained <- PeakSegFPOPLog(count.vec, penalty = pen)
      segments_con <- sum(constrained$mean.vec > 0, na.rm=TRUE)
      cost_con <- constrained$cost.mat[length(count.vec)*2-1]
      
      results <- rbind(results, data.table(
        penalty = pen,
        model = "Constrained",
        segments = segments_con,
        cost = cost_con
      ))
    }
  }
  
  return(results)
}

# Test with different data sizes
data_sizes <- c(100, 500, 1000)
all_results <- data.table()

for(size in data_sizes) {
  cat("Testing with data size:", size, "\n")
  
  # Generate data with 3 peaks
  data <- generate_synthetic_data(size, 3)
  
  # Define penalty range based on data size
  penalty_range <- 10^seq(0, 3, by=0.5)
  
  # Compare models
  results <- compare_models(data$count.vec, penalty_range)
  results$data_size <- size
  
  all_results <- rbind(all_results, results)
}

# Visualize results
p1 <- ggplot(all_results, aes(x=penalty, y=segments, color=model)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~data_size) +
  labs(title="Number of Segments vs Penalty", 
       x="Penalty (log scale)", 
       y="Number of Segments")

p2 <- ggplot(all_results, aes(x=penalty, y=cost, color=model)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~data_size) +
  labs(title="Model Cost vs Penalty", 
       x="Penalty (log scale)", 
       y="Cost")

print(p1)
print(p2)

# Save results for documentation
pdf("model-comparison-results.pdf", width=10, height=8)
print(p1)
print(p2)
dev.off()

# Print summary statistics
cat("\n=== Model Comparison Summary ===\n")
for(size in data_sizes) {
  cat("\nData size:", size, "\n")
  size_results <- all_results[data_size == size]
  
  # Compute average segments and costs
  avg_results <- size_results[, .(
    avg_segments = mean(segments),
    avg_cost = mean(cost)
  ), by=.(model)]
  
  print(avg_results)
  
  # Compare costs
  if("Constrained" %in% unique(size_results$model)) {
    cost_diff <- size_results[model == "Constrained", cost] - 
                 size_results[model == "Unconstrained", cost]
    cat("Average cost difference (Constrained - Unconstrained):", mean(cost_diff), "\n")
    cat("Percentage of cases where unconstrained has lower cost:", 
        mean(cost_diff > 0) * 100, "%\n")
  }
}
