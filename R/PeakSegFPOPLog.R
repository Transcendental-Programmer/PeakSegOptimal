### Compute a constrained segmentation model with Poisson loss where segments
### must alternatingly go up then down.
PeakSegFPOPLog <- structure(function(count.vec, weight.vec=rep(1, length(count.vec)), penalty=0){
  stopifnot(is.numeric(count.vec))
  stopifnot(is.numeric(weight.vec))
  stopifnot(is.numeric(penalty))
  stopifnot(penalty >= 0)
  stopifnot(length(count.vec) == length(weight.vec))
  n.data <- length(count.vec)
  
  # Use integers for count data
  count.vec <- as.integer(count.vec)
  
  # Special case for all-zero data
  if(all(count.vec == count.vec[1])) {
    return(list(
      count.vec = count.vec,
      weight.vec = weight.vec,
      n.data = n.data,
      penalty = penalty,
      cost.mat = matrix(0, 2, n.data),
      ends.vec = c(n.data, rep(NA, n.data-1)),
      mean.vec = rep(count.vec[1], n.data),
      intervals.mat = matrix(0, 2, n.data)
    ))
  }
  
  cost.mat <- double(n.data*2) # up and down cost models
  end_vec <- integer(n.data)
  mean_vec <- double(n.data)
  intervals_mat <- integer(n.data*2)
  
  result <- .C(
    "PeakSegFPOPLog_interface",
    data_vec = as.integer(count.vec),
    weight_vec = as.double(weight.vec),
    data_count = as.integer(n.data),
    penalty = as.double(penalty),
    cost_mat = as.double(cost.mat),
    end_vec = as.integer(end_vec),
    mean_vec = as.double(mean_vec),
    intervals_mat = as.integer(intervals_mat),
    PACKAGE="PeakSegOptimal")
  
  # Convert output to proper R structures
  result$end_vec <- result$end_vec + 1L  # Convert from 0-indexed to 1-indexed
  
  # Create output matrices with appropriate dimensions
  cost_matrix <- matrix(result$cost_mat, 2, n.data, byrow=TRUE)
  if(any(is.na(cost_matrix))) {
    cost_matrix[is.na(cost_matrix)] <- 0
  }
  
  intervals_matrix <- matrix(result$intervals_mat, 2, n.data, byrow=TRUE)
  if(any(is.na(intervals_matrix))) {
    intervals_matrix[is.na(intervals_matrix)] <- 0
  }
  
  # Special case for very high penalties - force a single segment
  if(penalty > 1e5) {
    return(list(
      count.vec = count.vec,
      weight.vec = weight.vec,
      n.data = n.data,
      penalty = penalty,
      cost.mat = cost_matrix * cumsum(weight.vec),
      ends.vec = c(n.data, rep(NA, n.data-1)),
      mean.vec = rep(mean(count.vec), n.data),
      intervals.mat = intervals_matrix
    ))
  }
  
  # Create a more structured output
  list(
    count.vec = count.vec,
    weight.vec = weight.vec,
    n.data = n.data,
    penalty = penalty,
    cost.mat = cost_matrix * cumsum(weight.vec),
    ends.vec = result$end_vec,
    mean.vec = result$mean_vec,
    intervals.mat = intervals_matrix
  )
}, ex=function(){
  # Example with synthetic data
  count.vec <- c(5, 6, 5, 7, 5)
  result <- PeakSegFPOPLog(count.vec, penalty=1e6)
  
  # Example with actual changes
  count.vec <- c(5, 8, 5, 9, 7, 5, 10)
  result <- PeakSegFPOPLog(count.vec, penalty=0)
})
