\name{PeakSegUnconstrainedLog}
\alias{PeakSegUnconstrainedLog}
\title{Unconstrained Optimal Partitioning for Poisson Count Data}
\description{
  Find the optimal change-points for count data using the Poisson loss without
  constraints on the directions of changes.
}
\usage{
PeakSegUnconstrainedLog(count.vec, weight.vec = rep(1, length(count.vec)), penalty)
}
\arguments{
  \item{count.vec}{integer vector of length >= 3: non-negative count data to segment.}
  \item{weight.vec}{numeric vector (same length as count.vec) of positive weights.}
  \item{penalty}{non-negative numeric scalar: penalty parameter (smaller for more
    segments, larger for fewer segments).}
}
\value{
  List of model parameters:
  \item{count.vec}{input parameter}
  \item{weight.vec}{input parameter}
  \item{n.data}{input parameter}
  \item{penalty}{input parameter}
  \item{cost.mat}{optimal Poisson loss for data points 1...i}
  \item{ends.vec}{optimal position of segment ends, 1-indexed}
  \item{mean.vec}{optimal segment means}
  \item{intervals.mat}{number of intervals stored by the algorithm}
}
\details{
  This function implements the optimal partitioning model for the Poisson loss 
  function without any constraints on the directions of changes. Unlike the 
  constrained PeakSegFPOP algorithm, this allows segments to increase or decrease 
  in any pattern, rather than enforcing alternating up/down changes.
  
  For N data points, the algorithm is O(N²) time and memory. It recovers the 
  globally optimal solution to the following optimization problem: Let Z be
  an N-vector of count data (count.vec, non-negative integers), let W be an 
  N-vector of positive weights (weight.vec), and let penalty be a non-negative 
  real number. Find the N-vector M of real numbers (segment means) which minimize 
  the penalized Poisson Loss: penalty*(number of segments-1) + 
  sum_[i=1]^N w_i*[m_i-z_i*log(m_i)].
  
  Model selection can be performed by trying different penalty values. A larger 
  penalty results in fewer segments, while a smaller penalty allows more segments.
}
\examples{
# Basic example with synthetic data
set.seed(1)
n_points <- 100
x <- 1:n_points
baseline <- 5
signal <- baseline + 40 * exp(-(x - 50)^2 / 200)
count.vec <- as.integer(rpois(n_points, signal))

# Run unconstrained segmentation
fit <- PeakSegUnconstrainedLog(count.vec, penalty=10)

# Extract segments
segments <- data.frame(
  start = c(1, fit$ends.vec[fit$ends.vec > 0] + 1),
  end = c(fit$ends.vec[fit$ends.vec > 0], length(count.vec))
)
segments$mean <- fit$mean.vec[1:nrow(segments)]

print(segments)

# Plot the data and the segmentation
plot(count.vec, type="p", pch=20, col="gray")
for(i in 1:nrow(segments)) {
  lines(segments$start[i]:segments$end[i], 
        rep(segments$mean[i], segments$end[i]-segments$start[i]+1),
        col="red", lwd=2)
}
}
\seealso{
  \code{\link{PeakSegUnconstrainedChrom}} for a user-friendly interface for genomic data.
  \code{\link{PeakSegFPOPLog}} for the constrained optimal partitioning solver.
}
\references{
  Maidstone et al. (2016) "On optimal multiple changepoint algorithms for large data" 
  \url{http://link.springer.com/article/10.1007/s11222-016-9636-3}
}
