# PeakSegOptimal

The **PeakSegOptimal** repository now includes several important modifications and new features that improve the range of segmentation models available for count data. In particular, we have implemented an unconstrained optimal partitioning solver for the Poisson loss model that does not enforce the up-down (alternating) constraints required by the previous FPOP solver. This new solver computes the optimal segmentation without restrictions on the direction of changes, making it more flexible when the data do not follow a strict peak–valley structure.

Below is an overview of the changes made, details about usage, commands, expected outputs, time complexities, and test instructions.

---

## Summary of Changes

- **New Unconstrained Solver:**
  - Implemented a new method `PiecewisePoissonLossLog::set_to_unconstrained_min_of` in the `funPieceListLog` module.
  - Added a new C++ function called `PeakSegUnconstrainedLog` which uses dynamic programming (DP) and cumulative sums to solve the optimal partitioning problem with Poisson loss _without any constraints_.
  - Unlike the constrained solver (implemented in `PeakSegFPOPLog.cpp`), which computes two sets (up and down) of cost functions (resulting in an N×2 array), the unconstrained solver produces a single state cost array of size N.
  - The unconstrained algorithm directly computes a constant piece representing the minimum cost function, leading to a final segmentation that uses dynamic programming with a double loop—yielding O(n²) worst-case time complexity.

- **Interface Updates:**
  - The `interface.cpp` file now includes an R interface for the unconstrained solver (`PeakSegUnconstrainedLog_interface`). This makes the new solver available in R with a call like `PeakSegUnconstrainedLog(...)`.
  - The new interface functions ensure that if the inputs do not allow a valid segmentation (e.g., all counts are equal to zero), an appropriate error is thrown.

- **Testing:**
  - Extensive `testthat` tests were added under the `tests/testthat` directory to compare the behavior of the unconstrained model against expected outcomes:
    - Basic functionality ensuring that without penalty every data point is segmented.
    - Edge cases, such as constant values, very high penalties (yielding a single segment), and weights affecting the segmentation.
    - A direct comparison with the constrained solver (`PeakSegFPOP`), ensuring that the unconstrained cost is lower or equal.
  - The tests also check that the unconstrained solver allows consecutive similar changes (up–up or down–down), which is not possible with the constrained model.

- **Documentation:**
  - This README summarizes the changes along with detailed usage instructions.
  - The implementation was inspired by previous work (see Maidstone et al. 2016) and the code in `PeakSegFPOPLog.cpp`.
  - For further details on the original segmentation approach, see the included documentation reference block (also appended as an Rd file).

---

## Unconstrained Segmentation Model

This implementation introduces an optimal partitioning algorithm for the Poisson loss without enforcing up–down constraints. It uses dynamic programming with efficient backtracking and has been benchmarked against the constrained FPOP method.

---

## How to Use

### In R

The new optimal partitioning solver for the unconstrained Poisson loss can be called directly using the R function:

```r
# Example usage for unconstrained segmentation:
count.vec <- as.integer(c(1, 10, 14, 13))
fit <- PeakSegUnconstrainedLog(count.vec, penalty = 0.1)

# Extract segmentation results
segments <- data.frame(
  start = c(1, fit$ends.vec[fit$ends.vec > 0] + 1),
  end = c(fit$ends.vec[fit$ends.vec > 0], length(count.vec))
)
segments$mean <- fit$mean.vec[1:nrow(segments)]
print(segments)
```

### R Interface and Commands

- **Constrained Solver:** Available via `PeakSegFPOPLog` (and its corresponding R interface function `PeakSegFPOP`).
- **Unconstrained Solver:** Available via `PeakSegUnconstrainedLog`.  
- The underlying C/C++ functions communicate with R through the `.C` interface, and the results are returned in a list containing:
  - `cost.mat`: Optimal Poisson loss computed for each end point.
  - `ends.vec`: Optimal segmentation end points (1-indexed in R).
  - `mean.vec`: Segment means for the optimal partition.
  - `intervals.mat`: Number of intervals stored by the algorithm.

All results (segmentation parameters and losses) are stored in these arrays and can be further processed or plotted.

### Additional Setup, Testing, Benchmarking & Model Selection Instructions

1. **Development Setup**
   - Build package:  
     ```r
     devtools::build()
     ```
   - Install package:  
     ```r
     devtools::install()
     ```
   - Check package:  
     ```r
     devtools::check()
     ```

2. **Running Tests**
   - Run tests using testthat:  
     ```r
     library(testthat)
     test_dir("tests/testthat")
     ```
   - Or run tests with devtools:  
     ```r
     devtools::test()
     ```

3. **Benchmarking & Model Selection**
   - Run the benchmarking script:  
     ```r
     source("tests/benchmark.R")
     ```
   - For model selection, if using a dedicated script (e.g., model-selection-compare.R):  
     ```r
     source("tests/model-selection-compare.R")
     ```
   - These scripts generate synthetic data, compare constrained vs. unconstrained segmentation, and create runtime and cost comparison plots.

4. **Using the New Algorithm**
   - Use the R interface for the new unconstrained segmentation as follows:
     ```r
     count.vec <- as.integer(c(1, 10, 14, 13))
     fit <- PeakSegUnconstrainedLog(count.vec, penalty = 0.1)
     segments <- data.frame(
       start = c(1, fit$ends.vec[fit$ends.vec > 0] + 1),
       end = c(fit$ends.vec[fit$ends.vec > 0], length(count.vec))
     )
     segments$mean <- fit$mean.vec[1:nrow(segments)]
     print(segments)
     ```

### Running Tests

The repository now uses the **testthat** package. To run the tests, use the command from the repository root in R:

```r
library(testthat)
test_dir("tests/testthat")
```

This will run the suite of tests ensuring that:
- Basic functionality and edge cases are handled correctly.
- The unconstrained solver provides a segmentation that is at least as good as the constrained one (in terms of cost).
- Segmentation outputs match expected behavior compared to the Segmentor3IsBack package for model selection.

### Testing Unconstrained Algorithm

To run tests for the new unconstrained algorithm only, execute:
```bash
Rscript -e "library(testthat); test_file('tests/testthat/test-PeakSegUnconstrainedLog.R')"
```

### Benchmarking & Model Selection

A benchmarking script `model-selection-compare.R` is provided which:
- Generates synthetic data with multiple peaks.
- Varies the penalty parameter to perform model selection.
- Uses **ggplot2** to visualize the number of segments and cost versus penalty (on a log scale).
- Saves the visualization results to a PDF (`model-comparison-results.pdf`).

There is also a `benchmark.R` test script under the tests folder to compare runtimes for different data sizes. Note:

- The constrained FPOP algorithm is O(n) in time.
- The newly implemented unconstrained solver runs in O(n²) time due to its double loop in the dynamic programming stage.

Given this, while the unconstrained solver is more flexible, it may be slower for very large datasets.

---

## Background and Reference

The original challenge was inspired by a need to implement an optimal partitioning method to compute the best K-segment model for a single penalty parameter without constraints. The unconstrained solver was modeled after the FPOP algorithm described in [Maidstone et al. (2016)](http://link.springer.com/article/10.1007/s11222-016-9636-3?wt_mc%3Dinternal.event.1.SEM.ArticleAuthorOnlineFirst). In contrast to the constrained approach (with up and down states and N×2 cost functions), the unconstrained solver computes just one set of cost functions (of size N) by calculating a constant piece representing the overall minimum.

The Segmentor3IsBack package implements the segment neighborhood model (best models in 1,...,K segments) for the Poisson loss with no constraints, however no implementation was originally available for the optimal partitioning model (best model for a single penalty parameter). This repository extends the work by implementing the missing solver and verifying it against the Segmentor3IsBack output for model selection.

---

## Conclusion

This update to **PeakSegOptimal** provides users with a more general segmentation tool through the unconstrained optimal partitioning solver. Users can now perform model selection over a range of penalty values, compare outputs between constrained (FPOP) and unconstrained models, and apply segmentation to various types of count data. Remember that while the unconstrained solver offers increased flexibility, it has a higher computational cost (O(n²)) compared to the O(n) constrained solver.

For any further questions or issues, please refer to the [BugReports](https://github.com/tdhock/PeakSegOptimal/issues) page.

Happy Segmenting!
