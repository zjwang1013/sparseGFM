# sparseGFM: Sparse Generalized Factor Models with Multiple Penalty Functions

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

The `sparseGFM` package provides methods for fitting sparse generalized factor models with various penalty functions. The package is designed to handle high-dimensional data and can adapt to weak factor scenarios, making it suitable for a wide range of applications in statistics and machine learning. 

The package is now available on [CRAN](https://cran.r-project.org/package=sparseGFM) under the name **sparseGFM**.  

## Installation

The package is now available on [CRAN](https://cran.r-project.org/package=sparseGFM) under the name **sparseGFM**.  
On GitHub, the development version is hosted under the repository name **sparseGFM**.

### From CRAN (stable version)

```r
install.packages("sparseGFM")

# Load the package
library(sparseGFM)
```

### From GitHub (development version)

You can also install the development version of sparseGFM from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install sparseGFM from GitHub
devtools::install_github("zjwang1013/sparseGFM")

# Load the package
library(sparseGFM)
```

## Key Features

- **Sparse Loading Matrix Estimation**: Efficiently estimates row-sparse loading matrices
- **Multiple Penalty Functions**: Supports lasso, adaptive group lasso (aglasso), and other penalties
- **Weak Factor Adaptation**: Automatically adapts to scenarios with weak factor structures
- **Cross-Validation**: Built-in cross-validation for optimal parameter selection
- **Determine the Number of Factors**: Automatic selection of the optimal number of factors using multiple information criteria

## Functions Overview

### Core Functions
- `sparseGFM()`: Main function for fitting sparse generalized factor models
- `cv.sparseGFM()`: Cross-validation for lambda selection
- `facnum.sparseGFM()`: Information criterion-based selection of factor number

### Supporting Functions
- `add_identifiability()`: Apply identifiability constraints to factor/loading matrices
- `evaluate_performance()`: Evaluate variable selection performance metrics

## Penalty Functions

The package implements 12 different penalty functions:

| Penalty | Description | Adaptive Version |
|---------|-------------|------------------|
| `lasso` | L1 penalty | `alasso` |
| `SCAD` | Smoothly Clipped Absolute Deviation | `agSCAD` |
| `MCP` | Minimax Concave Penalty | `agMCP` |
| `group`/`glasso` | Group Lasso | `agroup`/`aglasso` |
| `gSCAD` | Group SCAD | `agSCAD` |
| `gMCP` | Group MCP | `agMCP` |

## Missing Data

The package automatically handles missing values.

## Quick Start

### Load the Package

```r
library(sparseGFM)
set.seed(123)
```

### Data Generation with Sparse Loading Matrix

```r
# Parameters
n <- 200          # number of observations
p <- 200          # number of variables
a_param <- 0.9    # sparsity parameter
q <- 2            # number of factors

# Generate sparse structure
s <- ceiling(p^a_param)  # number of non-zero loadings

# Generate factor matrix
FF <- matrix(runif(n * q, min = -3, max = 3), nrow = n, ncol = q)

# Generate sparse loading matrix (row-sparse)
BB <- rbind(matrix(runif(s * q, min = 1, max = 2), nrow = s, ncol = q),
            matrix(0, nrow = (p - s), ncol = q))

# Generate intercepts
alpha_true <- runif(p, min = -1, max = 1)

# Add identifiability constraints
ident_res <- add_identifiability(FF, BB, alpha_true)
FF0 <- ident_res$H
BB0 <- ident_res$B
alpha0 <- ident_res$mu

# Generate data matrix
mat_para <- FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha0))
x <- matrix(nrow = n, ncol = p)
for (i in 1:n) {
  for (j in 1:p) {
    x[i, j] <- rnorm(1, mean = mat_para[i, j])
  }
}

# True variable selection indicator
index_true <- numeric(p)
if (s > 0 && s <= p) {
  index_true[1:s] <- 1
}
```

## Examples

### Example 1: Cross-Validation with Adaptive Group Lasso

```r
# Perform cross-validation to select optimal lambda
cv_result <- cv.sparseGFM(x,
                          type = "continuous",
                          q = 2,
                          penalty = "aglasso",
                          C = 5,
                          lambda_range = seq(0.1, 1, by = 0.1),
                          verbose = FALSE)

# Extract optimal model
optimal_model <- cv_result$optimal_model
BB_hat <- optimal_model$BB_hat
FF_hat <- optimal_model$FF_hat
alpha_hat <- optimal_model$alpha_hat

# Evaluate model performance
mat_para_hat <- FF_hat %*% t(BB_hat) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha_hat))
relative_error <- base::norm((mat_para_hat - mat_para), type = "F") / base::norm(mat_para, type = "F")

print(paste("Optimal lambda:", cv_result$optimal_lambda))
print(paste("Relative estimation error:", round(relative_error, 4)))

# Variable selection performance
index_pred <- rep(1, p)
index_pred[optimal_model$index] <- 0
performance <- evaluate_performance(index_true, index_pred)
print(performance)
```

### Example 2: Direct Model Fitting

```r
# Fit sparse GFM with fixed lambda
result <- sparseGFM(x, 
                    type = "continuous", 
                    q = 2,
                    penalty = "aglasso",
                    lambda = 0.1,
                    C = 5)

# Extract results
BB_direct <- result$BB_hat
FF_direct <- result$FF_hat
alpha_direct <- result$alpha_hat

# View selected variables
selected_vars <- setdiff(1:p, result$index)
print(paste("Number of selected variables:", length(selected_vars)))
print(paste("Number of zero loadings:", length(result$index)))

# Calculate space angles for evaluation
BB_vcc <- eval.space(BB_direct, BB0)[1]
BB_tcc <- eval.space(BB_direct, BB0)[2]
FF_vcc <- eval.space(FF_direct, FF0)[1] 
FF_tcc <- eval.space(FF_direct, FF0)[2]

print(paste("Loading matrix VCC:", round(BB_vcc, 4)))
print(paste("Loading matrix TCC:", round(BB_tcc, 4)))
```

### Example 3: Determine the Number of Factors

```r
# Select optimal number of factors using multiple criteria
facnum_result <- facnum.sparseGFM(x,
                                  type = "continuous",
                                  q_range = 1:5,
                                  penalty = "aglasso",
                                  lambda_range = c(0.1),
                                  sic_type = "sic1",
                                  C = 6,
                                  verbose = FALSE)

# Extract information criteria results
df_dd <- facnum_result$df_dd
df_as <- facnum_result$df_as

# Get optimal factor numbers from different criteria
optimal_q_sic1 <- facnum_result$optimal_q
optimal_q_sic2 <- which.min(facnum_result$sic2)
optimal_q_sic3 <- which.min(facnum_result$sic3)
optimal_q_sic4 <- which.min(facnum_result$sic4)

print("Optimal number of factors:")
print(paste("SIC1:", optimal_q_sic1))
print(paste("SIC2:", optimal_q_sic2))
print(paste("SIC3:", optimal_q_sic3))
print(paste("SIC4:", optimal_q_sic4))

# Plot information criteria (if plotting functions are available)
plot(1:5, facnum_result$sic1, type = "b", col = "red", 
     xlab = "Number of Factors", ylab = "Information Criterion",
     main = "Determine the Number of Factors")
lines(1:5, facnum_result$sic2, type = "b", col = "blue")
lines(1:5, facnum_result$sic3, type = "b", col = "green")
lines(1:5, facnum_result$sic4, type = "b", col = "purple")
legend("topright", c("SIC1", "SIC2", "SIC3", "SIC4"), 
       col = c("red", "blue", "green", "purple"), lty = 1)
```

## Parameters

### Main Parameters

- **`x`**: Data matrix (n × p)
- **`type`**: Data type ("continuous" for Gaussian data)
- **`q`**: Number of factors
- **`penalty`**: Penalty type ("lasso", "aglasso", etc.)
- **`lambda`**: Regularization parameter
- **`C`**: Constraint constant

### Cross-Validation Parameters

- **`lambda_range`**: Range of lambda values for cross-validation
- **`verbose`**: Whether to print progress information

### Determine the Number of Factors Parameters

- **`q_range`**: Range of factor numbers to consider
- **`sic_type`**: Type of information criterion for primary selection

## Methodological Details

The sparseGFM algorithm employs:

1. **Initialization**: Using the GFM package for initial estimates
2. **Alternating minimization**: Iteratively updating factors (F) and loadings (B)
3. **Sparsity induction**: Applying various penalty functions to achieve variable selection
4. **Identifiability**: Ensuring unique solutions through SVD-based constraints
5. **Convergence**: Monitoring objective function changes until convergence

## Model Characteristics

The sparseGFM package is specifically designed to handle:

1. **High-dimensional data** where the number of variables may exceed the number of observations
2. **Sparse loading structures** with row-wise sparsity patterns
3. **Weak factor scenarios** where factors have relatively small eigenvalues
4. **Flexible penalty structures** including adaptive penalties that can provide better variable selection

The adaptive group lasso (aglasso) penalty is particularly effective for row-sparse loading matrices, as it can select entire rows (variables) rather than individual elements.

## Dependencies

- R (≥ 3.5.0)
- GFM
- MASS
- irlba
- stats

## Bug Reports and Issues

Please report any bugs or issues on the [GitHub Issues page](https://github.com/zjwang1013/sparseGFM/issues). When reporting, include:

1. A clear description of the issue
2. Reproducible code example
3. Your session information (`sessionInfo()`)
4. Any relevant error messages

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss proposed modifications.

## References

For more detailed information about the methodology and theoretical properties, please refer to the associated research papers and documentation. The related papers are currently under review; we will update the information as soon as they are accepted and published.

## License

This package is licensed under GPL-3. See the [LICENSE](LICENSE) file for details.

## Issues and Support

For bug reports, feature requests, or questions, please visit the GitHub repository issues page.
