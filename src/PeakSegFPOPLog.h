int PeakSegFPOPLog
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 // the following matrices are for output.
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat);

int PeakSegUnconstrainedLog
(int *data_vec, double *weight_vec,
 int data_count, double penalty,
 // Output arrays: here we use a single cost function (size data_count x 1).
 double *cost_mat,
 int *end_vec,
 double *mean_vec,
 int *intervals_mat);
