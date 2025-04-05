library(PeakSegOptimal)
library(data.table)
library(ggplot2)

# Generate or load real-world test data
set.seed(123)
n_points <- 100
# Create synthetic data with multiple peaks
x <- 1:n_points
peak_positions <- c(20, 50, 80)
peak_heights <- c(30, 50, 40)
baseline <- 5
signal <- baseline
for(i in 1:length(peak_positions)) {
  signal <- signal + peak_heights[i] * exp(-(x - peak_positions[i])^2 / 50)
}
count.vec <- as.integer(rpois(n_points, signal))

# Run both algorithms for comparison
system.time(
  unconstrained_result <- PeakSegUnconstrainedLog(count.vec, penalty=15)
)

# For constrained model, use the R interface FPOP function
system.time(
  constrained_result <- PeakSegFPOP(count.vec, rep(1, length(count.vec)), penalty=15)
)

# Prepare visualization data
dt <- data.table(
  position = 1:n_points,
  count = count.vec
)

# Extract unconstrained results
unc_segments <- data.frame(
  start = c(1, unconstrained_result$ends.vec[unconstrained_result$ends.vec > 0] + 1),
  end = c(unconstrained_result$ends.vec[unconstrained_result$ends.vec > 0], length(count.vec))
)
unc_segments$mean <- unconstrained_result$mean.vec[1:nrow(unc_segments)]

# Extract constrained results 
con_segments <- data.frame(
  start = c(1, constrained_result$ends.vec[constrained_result$ends.vec > 0] + 1),
  end = c(constrained_result$ends.vec[constrained_result$ends.vec > 0], length(count.vec))
)
con_segments$mean <- constrained_result$mean.vec[1:nrow(con_segments)]

# Add segment means to data table
dt$unconstrained <- NA_real_
dt$constrained <- NA_real_

# Fix loop for unconstrained results
for(i in 1:nrow(unc_segments)) {
  segment_range <- unc_segments$start[i]:unc_segments$end[i]
  # Make sure we don't go out of bounds
  segment_range <- segment_range[segment_range <= nrow(dt)]
  if(length(segment_range) > 0) {
    dt[segment_range, unconstrained := unc_segments$mean[i]]
  }
}

# Fix loop for constrained results
for(i in 1:nrow(con_segments)) {
  segment_range <- con_segments$start[i]:con_segments$end[i]
  # Make sure we don't go out of bounds
  segment_range <- segment_range[segment_range <= nrow(dt)]
  if(length(segment_range) > 0) {
    dt[segment_range, constrained := con_segments$mean[i]]
  }
}

# Plot the results
p <- ggplot(dt) +
  geom_point(aes(x=position, y=count), color="black", alpha=0.5) +
  geom_line(aes(x=position, y=unconstrained, color="Unconstrained"), size=1) +
  geom_line(aes(x=position, y=constrained, color="Constrained"), size=1, linetype="dashed") +
  scale_color_manual(values=c("Unconstrained"="blue", "Constrained"="red")) +
  labs(title="Comparison of Constrained vs Unconstrained Peak Segmentation",
       y="Count", color="Method") +
  theme_minimal()

print(p)

# Calculate stats for comparison
cat("Constrained model info:\n")
cat("- Number of segments:", nrow(con_segments), "\n")

cat("\nUnconstrained model info:\n")
cat("- Number of segments:", nrow(unc_segments), "\n")

# Check execution time for different data sizes
sizes <- c(100, 500, 1000)
times <- data.frame(
  size = sizes,
  unconstrained = numeric(length(sizes)),
  constrained = numeric(length(sizes))
)

for(i in 1:length(sizes)) {
  n <- sizes[i]
  x <- 1:n
  signal <- baseline + 50 * exp(-(x - n/2)^2 / (n/10))
  counts <- as.integer(rpois(n, signal))
  
  t1 <- system.time(PeakSegUnconstrainedLog(counts, penalty=15))
  times$unconstrained[i] <- t1[3]
  
  t2 <- system.time(PeakSegFPOP(counts, rep(1, length(counts)), penalty=15))
  times$constrained[i] <- t2[3]
  
  cat("Size", n, "- Unconstrained:", t1[3], "s, Constrained:", t2[3], "s\n")
}

# Plot timing comparison
timing_plot <- ggplot(times, aes(x=size)) +
  geom_line(aes(y=unconstrained, color="Unconstrained")) +
  geom_line(aes(y=constrained, color="Constrained")) +
  labs(title="Runtime Comparison", x="Data Size", y="Execution Time (s)",
       color="Method") +
  theme_minimal()

print(timing_plot)
