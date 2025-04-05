/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceListLog.h"
#include <math.h>
#include <R.h> // for Rprintf

int IsotonicRegressionNormal
(double *data_vec, double *weight_vec, int data_count,
 double penalty,
 // output parameters:
 double *cost_vec, //data_count 
 int *end_vec,     //data_count
 double *mean_vec) //data_count
{
  // Initialize with min/max values for the data range
  double min_mean = -INFINITY, max_mean = INFINITY;
  for(int i=0; i<data_count; i++){
    if(data_vec[i] < min_mean || min_mean == -INFINITY){
      min_mean = data_vec[i];
    }
    if(max_mean < data_vec[i] || max_mean == INFINITY){
      max_mean = data_vec[i];
    }
  }
  
  // For numerical stability, if all data points are the same
  if(min_mean == max_mean){
    for(int i=0; i<data_count; i++){
      mean_vec[i] = min_mean;
      end_vec[i] = 0;
      cost_vec[i] = 0.0;
    }
    return 0;
  }
  
  // Buffer for min/max to avoid numerical issues
  double range = max_mean - min_mean;
  min_mean -= range * 0.1;
  max_mean += range * 0.1;
  
  // Create the cost model vector
  std::vector<PiecewiseNormalLoss> cost_model_vec(data_count);
  
  // Initialize first segment with the first data point
  PiecewiseNormalLoss *cost_model = &cost_model_vec[0];
  double weight = weight_vec[0];
  double observation = data_vec[0];
  
  // For normal loss, quadratic=weight, linear=-2*weight*observation, constant=weight*observation^2
  cost_model->piece_list.emplace_back(
    weight, -2*weight*observation, weight*observation*observation,
    min_mean, max_mean, -1, 0.0);
  
  // Fill in cost_vec, end_vec for the first data point
  double best_cost, best_mean;
  int best_end;
  double prev_mean;
  cost_model->Minimize(&best_cost, &best_mean, &best_end, &prev_mean);
  cost_vec[0] = best_cost;
  mean_vec[0] = best_mean;
  end_vec[0] = -1;  // No previous changepoint
  
  // DP recursion: optimal cost up to data point i
  for(int i=1; i<data_count; i++){
    // Create a new cost model for optimizing up to position i
    PiecewiseNormalLoss *prev_cost_model = &cost_model_vec[i-1];
    PiecewiseNormalLoss *new_cost_model = &cost_model_vec[i];
    
    // Normal cost terms for the current data point
    weight = weight_vec[i];
    observation = data_vec[i];
    PiecewiseNormalLoss data_cost;
    data_cost.piece_list.emplace_back(
      weight, -2*weight*observation, weight*observation*observation,
      min_mean, max_mean, -1, 0.0);
    
    // Case 1: No new segment - extend the previous model
    PiecewiseNormalLoss no_change_model = *prev_cost_model;
    no_change_model.add(
      weight, -2*weight*observation, weight*observation*observation);
    
    // Case 2: New segment starts here - enforce isotonic constraint
    PiecewiseNormalLoss min_prev_cost;
    // Use min_less_of to enforce isotonic constraint
    min_prev_cost.set_to_min_less_of(prev_cost_model, 0);
    min_prev_cost.set_prev_seg_end(i-1);
    
    // Apply penalty for the new changepoint
    if (penalty > 0) {
      min_prev_cost.add(0.0, 0.0, penalty);
    }
    
    // Add the data cost to the changepoint model
    PiecewiseNormalLoss change_model = min_prev_cost;
    change_model.add(
      weight, -2*weight*observation, weight*observation*observation);
    
    // Compare the no-change and change models to find minimum cost
    new_cost_model->piece_list.clear();
    double min_cost_nochange, min_mean_nochange;
    int end_nochange;
    double prev_mean_nochange;
    
    double min_cost_change, min_mean_change;
    int end_change;
    double prev_mean_change;
    
    no_change_model.Minimize(&min_cost_nochange, &min_mean_nochange, 
                           &end_nochange, &prev_mean_nochange);
    change_model.Minimize(&min_cost_change, &min_mean_change,
                        &end_change, &prev_mean_change);
    
    if(min_cost_change < min_cost_nochange){
      *new_cost_model = change_model;
      cost_vec[i] = min_cost_change;
      mean_vec[i] = min_mean_change;
      end_vec[i] = end_change;
    } else {
      *new_cost_model = no_change_model;
      cost_vec[i] = min_cost_nochange;
      mean_vec[i] = min_mean_nochange;
      end_vec[i] = end_nochange;
    }
  }
  
  return 0;
}
