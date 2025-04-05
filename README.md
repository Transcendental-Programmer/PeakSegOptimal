# PeakSegOptimal - Medium Test Branch

This branch implements rigorous testing for the unconstrained optimal partitioning algorithm with Poisson loss, focusing on challenging edge cases and numerical stability.

## Environment Setup

### Setting Up with renv

1. Clone the repository and checkout the master branch:
   ```bash
   git clone https://github.com/Transcendental-Programmer/PeakSegOptimal.git
   cd PeakSegOptimal
   git checkout master
   ```

2. Initialize the renv environment:
   ```r
   install.packages("renv")
   renv::init()
   ```

3. Install required packages:
   ```r
   renv::install("devtools")
   renv::install("testthat")
   renv::install("data.table")
   renv::install("ggplot2")
   ```

### Installing the Package

Build and install the package from source:

```r
devtools::build()
devtools::install()
```

Verify the installation was successful:

```r
library(PeakSegOptimal)
```

## Running Tests

Follow these steps to run the tests and examine the results:

1. **Run all tests**: This will execute all test files including the hard tests.
   ```r
   library(testthat)
   test_dir("tests/testthat")
   ```
   You should see test results showing passed, failed, and skipped tests.

2. **Run specific unconstrained segmentation tests**:
   ```r
   test_file("tests/testthat/test-PeakSegUnconstrainedLog.R")
   ```
   This will show detailed results for just the unconstrained tests.

3. **Test specific edge cases**: If you want to focus on particular edge cases:
   ```r
   test_file("tests/testthat/test-PeakSegUnconstrainedLog.R", filter="edge cases")
   ```
   This will only run tests with "edge cases" in their description.

4. **Benchmark performance**:
   ```r
   source("tests/benchmark.R")
   ```
   This will generate performance comparison plots and print metrics.

5. **Model selection comparison**:
   ```r
   source("model-selection-compare.R")
   ```
   This will compare model selection behavior across different penalty values.

## Output Files and Results

When running the tests and benchmark scripts, you'll find outputs in these locations:

1. **Test Results**: When using `test_dir()` or `test_file()`, test results are displayed in the R console. Detailed test logs can be found in:
   ```
   tests/testthat/test-results/
   ```

2. **Benchmark Results**: Running the benchmark script generates visualizations saved at:
   ```
   tests/benchmark-results/timing-comparison.pdf
   tests/benchmark-results/segmentation-comparison.pdf
   ```

3. **Model Selection Outputs**: The model selection script saves plots and data to:
   ```
   model-selection-results/model-comparison-plot.pdf
   model-selection-results/segment-counts.csv
   ```

4. **Function Examples**: Example outputs from the included example functions are displayed in the R console and can be saved to your working directory by assigning the plot outputs.

## Approach

The hard_test branch extends the unconstrained optimal partitioning algorithm with:

1. **Robust dynamic programming**: The implementation uses a more numerically stable approach to backtracking that handles cases where conventional dynamic programming might fail.

2. **Edge case handling**: Special cases that could cause numerical issues (like zero counts, constant values, or extreme penalties) are explicitly handled.

3. **Comprehensive testing**: Tests include challenging scenarios that stress the numerical stability of the segmentation algorithms.

4. **Performance validation**: Benchmarks ensure the unconstrained solver performs as expected even in difficult cases.

The core algorithm still uses dynamic programming with O(n²) time complexity, but with improved numerical stability.

## Implementation for Medium Test Problem

### Medium Test Challenge
The Segmentor3IsBack package implements the segment neighborhood model (best models in 1,...,K segments) for the Poisson loss and no constraints, but there is no implementation available for the optimal partitioning model (best K-segment model for a single penalty parameter, without computing the models with 1,...,K-1 segments). The goal of this test is to modify the code in the `PeakSegOptimal` package, in order to implement a solver for the optimal partitioning problem with Poisson loss and no constraints. Begin by studying [PeakSegFPOPLog.cpp](https://github.com/tdhock/coseg/blob/master/src/PeakSegFPOPLog.cpp) which implements the optimal partitioning model for the Poisson loss and the up-down constraints. There are two states in this model, up and down.  Since the up-down constrained solver has two states, there are N x 2 optimal cost functions to compute (`cost_model_mat` is of dimension `data_count*2`). The cost of being in the up state at `data_i` is `cost_model_mat[data_i]` and the cost of being in the down state is `cost_model_mat[data_i+data_count]`. The `min_prev_cost.set_to_min_less_of(down_cost_prev)` method enforces a non-decreasing constraint between adjacent segment means, for the state change down->up. Analogously, the `PiecewisePoissonLossLog::set_to_min_more_of` method enforces a non-increasing constraint for the state change up->down. To implement the unconstrained solver, you just need to implement a new `PiecewisePoissonLossLog::set_to_unconstrained_min_of` method that computes the minimum constant cost function (one `PoissonLossPieceLog` object), and then uses that to compute the N x 1 array of optimal cost functions (`cost_model_mat`). Read about the FPOP algorithm in [Maidstone et al 2016](http://link.springer.com/article/10.1007/s11222-016-9636-3?wt_mc%3Dinternal.event.1.SEM.ArticleAuthorOnlineFirst) for more info. When you are done with your implementation, check your work by comparing with the output of `Segmentor3IsBack::Segmentor(model=1)`. Perform model selection yourself for a range of penalty parameters. Using `testthat`, write some test cases which make sure your function gives the same exact model as the corresponding selected Segmentor model.

### File-by-File Solution Implementation

Our implementation solves the medium test problem across several files:

- **src/funPieceListLog.cpp/h**: 
  - Added the key method `set_to_unconstrained_min_of()` which computes the minimum constant cost function (one `PoissonLossPieceLog` object)
  - This is the core method that removes the constraints between segment means

- **src/PeakSegUnconstrainedLog.cpp**:
  - Implements the dynamic programming algorithm for the optimal partitioning problem with no constraints
  - Uses cumulative sums and backtracking to efficiently compute and decode segments
  - Handles special cases (high/low penalties, constant values) for numerical stability
  - Differs from constrained version by using a simpler O(n²) algorithm with single N x 1 cost array

- **src/interface.cpp**:
  - Added `PeakSegUnconstrainedLog_interface` to connect the C++ implementation with R
  - Provides error handling consistent with other functions in the package

- **R/PeakSegUnconstrainedLog.R**:
  - User-friendly R interface to the C++ implementation
  - Provides proper data validation and conversion between R and C++ data structures
  - Also includes a genomic data interface (`PeakSegUnconstrainedChrom`)

- **tests/testthat/test-PeakSegUnconstrainedLog.R**:
  - Tests basic functionality, edge cases, and weighted loss
  - Crucially compares unconstrained to constrained models to verify the unconstrained version produces equal or lower cost
  - Verifies that consecutive up-up or down-down changes are possible in the unconstrained model

- **model-selection-compare.R and tests/benchmark.R**:
  - Implements model selection across penalty values
  - Benchmarks performance comparing constrained vs. unconstrained algorithms
  - Generates visualizations of results for documentation and analysis

## File Changes

### Core Implementation

- **src/PeakSegUnconstrainedLog.cpp**: 
  - New implementation of the unconstrained optimal partitioning algorithm
  - Handles special cases for extreme penalty values and numerical stability
  - Implements efficient backtracking for segment reconstruction

- **src/funPieceListLog.cpp/h**:
  - Added `set_to_unconstrained_min_of()` method for computing constant pieces
  - Enhanced root-finding and numerical stability

- **src/interface.cpp**:
  - Added `PeakSegUnconstrainedLog_interface` for R integration
  - Error handling for edge cases

### R Interface

- **R/PeakSegUnconstrainedLog.R**:
  - User-facing functions for unconstrained segmentation
  - Implementation of `PeakSegUnconstrainedLog()` and `PeakSegUnconstrainedChrom()`
  - Proper conversion between C++ and R data structures

### Tests and Benchmarks

- **tests/testthat/test-PeakSegUnconstrainedLog.R**:
  - Comprehensive edge case testing
  - Comparison with constrained models
  - Validation of mathematical properties

- **tests/benchmark.R**:
  - Performance comparison between constrained and unconstrained algorithms
  - Tests across different data sizes

- **model-selection-compare.R**:
  - Assessment of model selection behavior with different penalty values

### Documentation

- **man/PeakSegUnconstrainedLog.Rd**:
  - Complete documentation of the new functions
  - Examples showing typical usage

## Contributing

For contributing to this branch:

1. Run the test suite to ensure stability: `devtools::test()`
2. Add new edge cases to the relevant test files
3. Document any numerical issues encountered

## License

This project is licensed under GPL-3 - see the LICENSE file for details.

## Acknowledgements

Based on the algorithms described in Maidstone et al. (2016) "On optimal multiple changepoint algorithms for large data".
