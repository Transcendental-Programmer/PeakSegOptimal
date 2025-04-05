/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <math.h>
#include <R.h> // for Rprintf
#include <algorithm>   // for std::reverse
#include "funPieceListLog.h"  // for ERROR_MIN_MAX_SAME

int PeakSegUnconstrainedLog
(int *data_vec, double *weight_vec, int data_count,
 double penalty,
 double *cost_mat, // data_count x 1.
 int *end_vec,     // data_count
 double *mean_vec, // data_count
 int *intervals_mat){ // data_count x 1.

  // Early check: if all counts are zero, error.
  int total = 0;
  for(int i=0; i<data_count; i++){
    total += data_vec[i];
  }
  if(total == 0){
    return ERROR_MIN_MAX_SAME;
  }

  // Initialize ALL output arrays to empty/default values
  for (int i = 0; i < data_count; i++){
    cost_mat[i] = NAN;
    end_vec[i] = -2;
    mean_vec[i] = NAN;
    intervals_mat[i] = 0;
  }

  // SPECIAL CASE: With huge penalty, just return a single segment with the overall mean
  if(penalty > 1e5) {
    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    for(int i = 0; i < data_count; i++){
      weighted_sum += data_vec[i] * weight_vec[i];
      weight_sum += weight_vec[i];
    }
    double overall_mean = weighted_sum / weight_sum;
    
    // Compute the Poisson loss for this single segment
    double loss = 0.0;
    for(int i = 0; i < data_count; i++){
      if(data_vec[i] > 0){
        loss += weight_vec[i] * (overall_mean - data_vec[i] * log(overall_mean));
      } else {
        loss += weight_vec[i] * overall_mean;
      }
    }
    
    // Store result: Only the first index gets the overall_mean; others are set to -1.
    cost_mat[data_count - 1] = loss;
    mean_vec[0] = overall_mean;
    end_vec[0] = data_count - 1;
    intervals_mat[0] = 1;
    for (int i = 1; i < data_count; i++){
      mean_vec[i] = -1;
      end_vec[i] = -2;
      intervals_mat[i] = 0;
      cost_mat[i] = NAN; // already initialized, but ensure unused entries remain dummy.
    }
    return 0;
  }

  // Compute cumulative sums for weights and weighted counts.
  std::vector<double> cum_weight(data_count,0.0);
  std::vector<double> cum_data(data_count,0.0);
  cum_weight[0] = weight_vec[0];
  cum_data[0] = data_vec[0]*weight_vec[0];
  for (int i=1; i<data_count; i++){
    cum_weight[i] = cum_weight[i-1] + weight_vec[i];
    cum_data[i] = cum_data[i-1] + data_vec[i]*weight_vec[i];
  }

  // DP arrays: dp[] and break_point[]
  std::vector<double> dp(data_count, INFINITY);
  std::vector<int> break_point(data_count, -1);
  
  // SPECIAL CASE: With near-zero penalty, create segments at every data point
  if(penalty < 1e-9){
    dp[0] = (data_vec[0] > 0) ? weight_vec[0]*(data_vec[0] - data_vec[0]*log((double)data_vec[0])) : 0.0;
    break_point[0] = 0;
    for(int i = 1; i < data_count; i++){
      double m = (data_vec[i] > 0) ? data_vec[i] : 1;
      double single_cost = (data_vec[i] > 0) ? weight_vec[i]*(m - m*log(m)) : 0.0;
      dp[i] = dp[i-1] + single_cost;
      break_point[i] = i;
    }
  } else {
    // Standard DP recurrence
    for (int i = 0; i < data_count; i++){
      for (int j = 0; j <= i; j++){
        double W = (j == 0) ? cum_weight[i] : (cum_weight[i] - cum_weight[j-1]);
        double Y = (j == 0) ? cum_data[i]   : (cum_data[i] - cum_data[j-1]);
        double m = (W > 0) ? (Y / W) : 0;
        double seg_cost = m * W - Y * log(m);
        double total_cost = seg_cost + (j > 0 ? dp[j-1] + penalty : 0);
        if(total_cost < dp[i]){
          dp[i] = total_cost;
          break_point[i] = j;
        }
      }
    }
  }
  
  // Decode segmentation using break_point
  std::vector<double> seg_means;
  std::vector<int> seg_ends;
  
  // Backtrack to get segments
  int i = data_count - 1;
  while(i >= 0){
    int j = break_point[i];
    double W = (j == 0) ? cum_weight[i] : (cum_weight[i] - cum_weight[j-1]);
    double Y = (j == 0) ? cum_data[i]   : (cum_data[i] - cum_data[j-1]);
    double m = (W > 0) ? (Y / W) : 0;
    seg_means.push_back(m);
    seg_ends.push_back(i);
    i = j - 1;
  }
  std::reverse(seg_means.begin(), seg_means.end());
  std::reverse(seg_ends.begin(), seg_ends.end());
  
  // Fill output arrays
  for(int i = 0; i < data_count; i++){
    cost_mat[i] = dp[i];
  }
  
  int num_seg = seg_means.size();
  for(int i = 0; i < num_seg; i++){
    mean_vec[i] = seg_means[i];
    end_vec[i] = seg_ends[i];
    intervals_mat[i] = 1;
  }
  
  return 0;
}