#ifndef ISOTONIC_REGRESSION_NORMAL_H
#define ISOTONIC_REGRESSION_NORMAL_H

int IsotonicRegressionNormal
(double *data_vec, double *weight_vec, int data_count,
 double penalty,
 // output parameters:
 double *cost_vec, //data_count 
 int *end_vec,     //data_count
 double *mean_vec); //data_count

#endif
