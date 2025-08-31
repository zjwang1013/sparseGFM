# sparseGFM: Sparse Generalized Factor Models with Multiple Penalty Functions

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

The **sparseGFM** package implements sparse generalized factor models for dimension reduction and variable selection in high-dimensional data. It provides a comprehensive framework for analyzing continuous, count, and binary data through appropriate generalized linear model frameworks with various sparsity-inducing penalties.

The package features automatic hyperparameter selection through cross-validation and information criteria, making it suitable for applications in genomics, economics, finance, and social sciences where interpretable dimension reduction is crucial.

## Key Features

- **Multiple data types**: Supports continuous (Gaussian), count (Poisson), and binary (Binomial) data
- **Comprehensive penalty functions**: Implements 12 different penalties including:
  - Standard penalties: Lasso, SCAD, MCP
  - Adaptive versions: Adaptive Lasso, Adaptive SCAD, Adaptive MCP
  - Group penalties: Group Lasso, Group SCAD, Group MCP
  - Adaptive group penalties: Adaptive Group Lasso, Adaptive Group SCAD, Adaptive Group MCP
- **Automatic parameter selection**: 
  - Cross-validation with BIC for regularization parameter (λ) selection
  - Information criteria for determining the number of factors (q)
- **Missing data support**: Automatic imputation for missing values
- **Efficient implementation**: Alternating minimization algorithm with identifiability constraints

## Installation

### From GitHub

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install sparseGFM from GitHub
devtools::install_github("zjwang1013/sparseGFM")

# Load the package
library(sparseGFM)
```

## Basic Usage

### Main Function: `sparseGFM()`

The primary function for fitting sparse generalized factor models:

```r
sparseGFM(x, type = "continuous", q = 2, penalty = "lasso", 
         lambda = 1, gam = 1, tau = NULL, mat_sd = 1, 
         delta = 1e-4, maxiter = 30, C = 5, verbose = TRUE)
```

**Parameters:**
- `x`: Data matrix (n × p)
- `type`: Data type - "continuous", "count", or "binary"
- `q`: Number of latent factors (default: 2)
- `penalty`: Penalty type (12 options available)
- `lambda`: Regularization parameter
- `gam`: Adaptive weight parameter for adaptive penalties
- `tau`: Shape parameter for SCAD/MCP penalties
- `C`: Constraint bound for stability (default: 5)
- `delta`: Convergence tolerance (default: 1e-4)
- `maxiter`: Maximum iterations (default: 30)

**Returns:** A list containing:
- `FF_hat`: Estimated factor matrix (n × q)
- `BB_hat`: Estimated loading matrix (p × q)
- `alpha_hat`: Estimated intercepts
- `obj_loglik`: Log-likelihood value
- `obj_pen`: Penalized objective value
- `index`: Indices of variables with zero loadings

## Examples

### Example 1: Basic Usage with Continuous Data

```r
library(sparseGFM)
set.seed(123)

# Generate continuous data
n <- 100  # number of observations
p <- 50   # number of variables
x <- matrix(rnorm(n * p), n, p)

# Fit sparse GFM with lasso penalty
result <- sparseGFM(x, type = "continuous", q = 3, 
                    penalty = "lasso", lambda = 0.5)

# View selected variables
selected_vars <- setdiff(1:p, result$index)
print(paste("Number of selected variables:", length(selected_vars)))
```

### Example 2: Cross-Validation for Lambda Selection

```r
# Automatic lambda selection via cross-validation
cv_result <- cv.sparseGFM(x, type = "continuous", q = 3,
                          penalty = "SCAD",
                          lambda_range = seq(0.1, 1, by = 0.1))

# Optimal lambda and model
print(paste("Optimal lambda:", cv_result$optimal_lambda))
optimal_model <- cv_result$optimal_model
```

### Example 3: Automatic Factor Number Selection

```r
# Determine optimal number of factors
facnum_result <- facnum.sparseGFM(x, type = "continuous",
                                  q_range = 1:6,
                                  penalty = "gSCAD",
                                  lambda_range = seq(0.1, 0.5, by = 0.05))

print(paste("Optimal number of factors:", facnum_result$optimal_q))
```

### Example 4: Count Data with Group Penalties

```r
# Generate count data
n <- 80; p <- 60
count_data <- matrix(rpois(n * p, lambda = 2), n, p)

# Fit with group lasso penalty
result <- sparseGFM(count_data, type = "count", q = 2,
                    penalty = "glasso", lambda = 0.3)

# Check sparsity pattern
print(table(rowSums(result$BB_hat^2) > 1e-6))
```

### Example 5: Binary Data with Adaptive Penalties

```r
# Generate binary data
binary_data <- matrix(rbinom(n * p, 1, 0.5), n, p)

# Fit with adaptive MCP penalty
result <- sparseGFM(binary_data, type = "binary", q = 2,
                    penalty = "agMCP", lambda = 0.2, gam = 0.5)
```

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

The package automatically handles missing values through mean imputation:
- Continuous data: Column means
- Count data: Rounded column means
- Binary data: Rounded column means (0 or 1)

Rows or columns with all missing values will trigger an error and should be removed before analysis.

## Methodological Details

The sparseGFM algorithm employs:
1. **Initialization**: Using the GFM package for initial estimates
2. **Alternating minimization**: Iteratively updating factors (F) and loadings (B)
3. **Sparsity induction**: Applying various penalty functions to achieve variable selection
4. **Identifiability**: Ensuring unique solutions through SVD-based constraints
5. **Convergence**: Monitoring Frobenius norm changes until convergence

## Citation

If you use sparseGFM in your research, please cite:

```
Wang, Z. (2024). sparseGFM: Sparse Generalized Factor Models with Multiple Penalty Functions. 
R package version 0.1.0. https://github.com/zjwang1013/sparseGFM
```

## Bug Reports and Issues

Please report any bugs or issues on the [GitHub Issues page](https://github.com/zjwang1013/sparseGFM/issues). When reporting, include:

1. A clear description of the issue
2. Reproducible code example
3. Your session information (`sessionInfo()`)
4. Any relevant error messages

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss proposed modifications.

## Author

**Zhijing Wang**  
Email: wangzhijing@sjtu.edu.cn  
Shanghai Jiao Tong University

## License

This package is licensed under GPL-3. See the [LICENSE](LICENSE) file for details.

## Dependencies

- R (≥ 3.5.0)
- GFM
- MASS
- irlba
- stats
