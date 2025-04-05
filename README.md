# PeakSegOptimal

This package implements optimal segmentation algorithms subject to constraints, including the new isotonic regression with normal loss.

## Environment Setup

### Prerequisites

- R (version 4.0.0 or later)
- C++ compiler (compatible with R)

### Installation

#### Using renv (recommended)

This project uses `renv` for reproducible environment management:

```r
# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv and restore dependencies
renv::init()
renv::restore()
```

#### Manual Installation

If you prefer not to use renv, you can install the required packages manually:

```r
install.packages(c("devtools", "testthat", "penaltyLearning", "ggplot2", "data.table"))
```

#### Installing the Development Version

To install the development version from this branch:

```r
devtools::install_github("tdhock/PeakSegOptimal", ref = "hard_test")
```

Or from your local clone:

```r
devtools::install()
```

## Running Tests

### Testing the Isotonic Regression Implementation

Run the complete test suite:

```r
devtools::test()
```

Run only the isotonic regression tests:

```r
devtools::test(filter = "isotonic")
```

You should see output similar to:

```
> devtools::test(filter = "isotonic-normal")
ℹ Testing PeakSegOptimal
✔ | F W  S  OK | Context
✔ |         14 | Isotonic Regression with Normal Loss                                       

══ Results ═════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 14 ]
> 
```

### Testing PeakSegFPOPLog

To test the PeakSegFPOPLog functionality:

```r
devtools::test(filter = "FPOPLog")
```

Expected output:

```
> devtools::test(filter = "FPOPLog")
ℹ Testing PeakSegOptimal
✔ | F W  S  OK | Context
✔ |          3 | FPOPLog

══ Results ═════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ]
```

### Manual Testing

You can manually test the isotonic regression function:

```r
library(PeakSegOptimal)

# Generate test data
set.seed(123)
data.vec <- rnorm(20)
data.vec[5:10] <- data.vec[5:10] + 2
data.vec[15:20] <- data.vec[15:20] + 3

# Run isotonic regression with different penalty values
result0 <- IsotonicRegression(data.vec, penalty=0) # No penalty
result5 <- IsotonicRegression(data.vec, penalty=5) # Medium penalty
result100 <- IsotonicRegression(data.vec, penalty=100) # High penalty

# Compare number of segments
nrow(result0)   # Many segments (no penalty)
nrow(result5)   # Fewer segments
nrow(result100) # Single segment

# Visualize results
plot(data.vec, type="p", pch=20)
lines(rep(result0$mean, result0$end - result0$start + 1), col="blue")
lines(rep(result5$mean, result5$end - result5$start + 1), col="red")
lines(rep(result100$mean, result100$end - result100$start + 1), col="green")
legend("topleft", legend=c("No penalty", "Medium penalty", "High penalty"), 
       col=c("blue", "red", "green"), lty=1)
```

## Project Challenge: Hard Test Solution

This project implements a regularized isotonic regression solver for normal loss, addressing the following challenge:

> ### Hard
> There is not yet an regularized isotonic regression solver for the normal loss ([issue](https://github.com/tdhock/coseg/issues/3)), and your goal in this test is to implement one. Like the unconstrained model described the Medium test, the regularized isotonic regression model also has only one state. So you can start by modifying the Medium test code, which should have a `cost_model_mat` which is N x 1. However the isotonic regression constraint means that all changes are non-decreasing, so you should use `set_to_min_less_of` instead of `set_to_unconstrained_min_of`. Now the difficult part: the existing code in the `coseg` package implements the Poisson loss via `class PoissonLossPieceLog`, but you need to implement another class for the Normal loss, `NormalLossPiece`. This class will need to declare different coefficients `Constant`, `Linear`, `Quadratic` for a function f(mean)=Constant + Linear*mean + Quadratic*mean^2. You will need to provide implementations for `get_smaller_root` and `get_larger_root` by using the `sqrt` function in `#include <math.h>`. It will be judged even better if you can get `PoissonLossPieceLog` and `NormalLossPiece` to inherit from the same base class with shared methods (that is the approach that the Segmentor package uses to implement several loss functions, and that is the approach that will be recommended for this GSOC project). Check your work by writing a `testthat` unit test to make sure that the model returned by your function with penalty=0 is the same as the model returned by the `isoreg` function (PAVA algorithm). Write another test that checks that the output model is the same as `Fpop` (when all changes are non-decreasing).

### Solution Components

Each file in the implementation addresses specific aspects of the challenge:

| File | Purpose | Challenge Component |
|------|---------|---------------------|
| **funPieceListLog.h** | Header file with class definitions | Defines the `LossPieceBase` base class that both `PoissonLossPieceLog` and `NormalLossPiece` inherit from; declares the `NormalLossPiece` class with `Quadratic`, `Linear`, and `Constant` coefficients |
| **funPieceListLog.cpp** | Implementation of loss functions | Implements `NormalLossPiece` methods including `get_smaller_root` and `get_larger_root` using `sqrt`; implements `set_to_min_less_of` for enforcing isotonic constraints |
| **IsotonicRegressionNormal.h** | Header for isotonic regression | Declares the function prototype for the solver |
| **IsotonicRegressionNormal.cpp** | C++ implementation | Implements the dynamic programming algorithm using `set_to_min_less_of` for isotonic constraints; handles penalty parameter for regularization |
| **IsotonicRegression.R** | R interface | Provides user-friendly access to the C++ solver; handles special cases and formatting of results |
| **NormalLoss.R** | Normal loss in R | Implements the weighted squared error loss function for testing and comparison |
| **interface.cpp** | C/R interface | Adds the `IsotonicRegressionNormal_interface` function to connect R with the C++ implementation |
| **test-isotonic-normal.R** | Comprehensive tests | Tests that the model matches `isoreg` (PAVA algorithm) when penalty=0; tests that the model matches FPOP with non-decreasing constraints |

## Approach

The implementation of isotonic regression with normal loss follows these key principles:

1. **Dynamic Programming**: The algorithm uses a dynamic programming approach to find the optimal segmentation under monotonicity constraints.

2. **Normal Loss Function**: Unlike the existing Poisson loss for count data, the normal loss is suitable for continuous data, implementing a squared error loss function.

3. **Regularization**: A penalty parameter allows controlling the number of segments, with higher penalties resulting in fewer segments.

4. **C++ Core with R Interface**: The core algorithm is implemented in C++ for efficiency, with a user-friendly R interface.

5. **Integration with Existing Framework**: The implementation follows the same pattern as other segmentation methods in the package, with consistent inputs and outputs.

## Implementation Details

### Added/Modified Files

#### R Files:
- **IsotonicRegression.R**: Main R interface for isotonic regression with normal loss
- **NormalLoss.R**: Implementation of the normal (squared error) loss function
- **PeakSegFPOPLog.R**: Fixed implementation for the PeakSegFPOPLog interface

#### C++ Files:
- **IsotonicRegressionNormal.cpp**: Core implementation of the isotonic regression algorithm
- **IsotonicRegressionNormal.h**: Header file declaring the function prototype
- **funPieceListLog.cpp**: Extended to include normal loss functionality
- **funPieceListLog.h**: Added base classes and normal loss pieces definitions
- **interface.cpp**: Added interface to the IsotonicRegressionNormal function

#### Test Files:
- **test-isotonic-normal.R**: Comprehensive tests for the isotonic regression implementation
- **test-FPOPLog.R**: Tests for PeakSegFPOPLog functionality

### Key Changes

1. **Base Class Architecture**:
   - Added `LossPieceBase` as a base class for different loss functions
   - Modified `PoissonLossPieceLog` to inherit from the base class
   - Added `NormalLossPiece` class for normal loss

2. **Normal Loss Implementation**:
   - Implemented quadratic loss functions
   - Added methods for finding roots, minimums, and computing costs
   
3. **Isotonic Constraint Handling**:
   - Implemented `set_to_min_less_of` for enforcing monotonicity constraints
   - Created specialized tests to ensure monotonicity

4. **Regularization**:
   - Added penalty parameter to control the number of change points
   - Implemented special cases for high penalties

5. **Integration with R**:
   - Created a comprehensive R interface
   - Added documentation and examples
   - Ensured compatibility with existing package functions

## License

This project is licensed under GPL-3.

## Contributing

Contributions are welcome! Please submit issues and pull requests through GitHub.
