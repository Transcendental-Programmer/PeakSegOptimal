PeakSegUnconstrainedLog <- structure(function
### Find the optimal change-points using the Poisson loss without
### constraints on the directions of changes. For N data points, the
### algorithm is O(N) time and memory. It recovers the globally
### optimal solution to the following optimization problem: Let Z be
### an N-vector of count data (count.vec, non-negative integers), let
### W be an N-vector of positive weights (weight.vec), and let penalty
### be a non-negative real number. Find the N-vector M of real numbers
### (segment means) which minimize the penalized Poisson Loss,
### penalty*(number of segments-1) + sum_[i=1]^N
### w_i*[m_i-z_i*log(m_i)].
(count.vec,
### integer vector of length >= 3: non-negative count data to segment.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL
### non-negative numeric scalar: penalty parameter (smaller for more
### segments, larger for fewer segments).
){
  n.data <- length(count.vec)
  stopifnot(3 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data==length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  
  # Initialize output arrays
  cost.mat <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.mat <- integer(n.data)
  
  result.list <- .C(
    "PeakSegUnconstrainedLog_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    cost.mat=as.double(cost.mat),
    ends.vec=as.integer(ends.vec),
    mean.vec=as.double(mean.vec),
    intervals.mat=as.integer(intervals.mat),
    PACKAGE="PeakSegOptimal")
  
  # Convert to 1-indexed segment ends for R
  result.list$ends.vec <- result.list$ends.vec+1L
  
  # Apply cumulative weights to costs
  result.list$cost.mat <- result.list$cost.mat * cumsum(weight.vec)
  
  result.list
### List of model parameters. count.vec, weight.vec, n.data, penalty
### (input parameters), cost.mat (optimal Poisson loss), ends.vec
### (optimal position of segment ends, 1-indexed), mean.vec (optimal
### segment means), intervals.mat (number of intervals stored by the
### algorithm).
}, ex=function(){
  
  # Basic example with synthetic data
  set.seed(1)
  n_points <- 100
  x <- 1:n_points
  baseline <- 5
  signal <- baseline + 40 * exp(-(x - 50)^2 / 200)
  count.vec <- as.integer(rpois(n_points, signal))
  
  # Run unconstrained segmentation
  fit <- PeakSegUnconstrainedLog(count.vec, penalty=10)
  
  # Examine the segmentation results
  segments <- data.frame(
    start = c(1, fit$ends.vec[fit$ends.vec > 0] + 1),
    end = c(fit$ends.vec[fit$ends.vec > 0], length(count.vec))
  )
  segments$mean <- fit$mean.vec[1:nrow(segments)]
  
  print(segments)
  
  # Plot the data and the segmentation
  library(ggplot2)
  plot_data <- data.frame(
    position = 1:length(count.vec),
    count = count.vec
  )
  
  # Create segmentation lines
  seg_lines <- data.frame()
  for(i in 1:nrow(segments)) {
    seg_lines <- rbind(seg_lines, data.frame(
      position = segments$start[i]:segments$end[i],
      mean = segments$mean[i]
    ))
  }
  
  ggplot(plot_data, aes(x=position, y=count)) +
    geom_point() +
    geom_line(data=seg_lines, aes(x=position, y=mean), color="red") +
    labs(title="Unconstrained Peak Segmentation", x="Position", y="Count") +
    theme_minimal()
})

PeakSegUnconstrainedChrom <- structure(function
### Find the optimal change-points using the Poisson loss and no constraints.
### This function is a user-friendly interface to the PeakSegUnconstrainedLog function.
(count.df,
### data.frame with columns count, chromStart, chromEnd.
 penalty=NULL
### non-negative numeric scalar: penalty parameter (smaller for more
### segments, larger for fewer segments).
){
  stopifnot(is.data.frame(count.df))
  n.data <- nrow(count.df)
  stopifnot(3 <= n.data)
  stopifnot(is.integer(count.df$chromStart))
  stopifnot(is.integer(count.df$chromEnd))
  stopifnot(is.integer(count.df$count))
  stopifnot(count.df$chromStart < count.df$chromEnd)
  stopifnot(0 <= count.df$chromStart)
  stopifnot(0 <= count.df$count)
  weight.vec <- with(count.df, chromEnd - chromStart)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  fit <- PeakSegUnconstrainedLog(count.df$count, weight.vec, penalty)
  
  # Extract segments
  is_change <- fit$ends.vec > 0
  break_vec <- fit$ends.vec[is_change]
  first <- c(1, break_vec+1)
  last <- c(break_vec, nrow(count.df))
  mean_vec <- fit$mean.vec[1:length(first)]
  
  list(
    segments=data.frame(
      mean=mean_vec,
      first,
      last,
      chromStart=count.df$chromStart[first],
      chromEnd=count.df$chromEnd[last],
      segments=length(first)),
    loss=data.frame(
      segments=length(first),
      penalized.loss=fit$cost.mat[n.data],
      feasible=TRUE  # Always feasible since there are no constraints
    )
  )
### List of data.frames: segments can be used for plotting the
### segmentation model, loss summarizes the penalized PoissonLoss.
}, ex=function(){
  # Create sample genomic data
  n <- 100
  chrom_df <- data.frame(
    count = as.integer(c(rpois(n/3, 5), rpois(n/3, 20), rpois(n/3, 5))),
    chromStart = seq(0, n-1) * 100,
    chromEnd = seq(1, n) * 100
  )
  
  # Run unconstrained segmentation
  fit <- PeakSegUnconstrainedChrom(chrom_df, 100)
  
  # Plot the results
  library(ggplot2)
  ggplot() +
    geom_step(aes(x=chromStart/1e3, y=count), data=chrom_df, color="grey") +
    geom_segment(aes(x=chromStart/1e3, y=mean, xend=chromEnd/1e3, yend=mean),
                 color="blue", data=fit$segments) +
    labs(title="Unconstrained Segmentation", 
         x="Position (kb)", 
         y="Coverage") +
    theme_minimal()
})
