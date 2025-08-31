#' @import GFM
#' @importFrom stats binomial cancor coef dbinom dnorm dpois gaussian glm glm.fit plogis poisson
#' @importFrom MASS ginv
#' @importFrom irlba irlba
NULL

#' Sparse Generalized Factor Model
#'
#' @description
#' Implements sparse generalized factor models with various penalty functions for dimension reduction
#' and variable selection in high-dimensional data with different data types (continuous, count, binary).
#'
#' @param x A numeric matrix of observations (n x p), where n is the number of observations
#'        and p is the number of variables
#' @param type Character string specifying the data type. Options are:
#'        \itemize{
#'          \item "continuous": Gaussian family for continuous data
#'          \item "count": Poisson family for count data
#'          \item "binary": Binomial family for binary data
#'        }
#' @param q Integer specifying the number of latent factors (default = 2)
#' @param penalty Character string specifying the penalty type for sparsity. Options include:
#'        \itemize{
#'          \item "lasso": L1 penalty
#'          \item "alasso": Adaptive lasso penalty
#'          \item "SCAD": Smoothly clipped absolute deviation penalty
#'          \item "MCP": Minimax concave penalty
#'          \item "group"/"glasso": Group lasso penalty
#'          \item "agroup"/"aglasso": Adaptive group lasso penalty
#'          \item "gSCAD": Group SCAD penalty
#'          \item "agSCAD": Adaptive group SCAD penalty
#'          \item "gMCP": Group MCP penalty
#'          \item "agMCP": Adaptive group MCP penalty
#'        }
#' @param lambda Numeric value for the penalty tuning parameter (default = 1)
#' @param gam Numeric value for the adaptive weight parameter in adaptive penalties (default = 1)
#' @param tau Numeric value for the shape parameter in SCAD/MCP penalties.
#'        Default is 3.7 for SCAD and 3 for MCP if not specified
#' @param mat_sd Standard deviation for continuous data (default = 1)
#' @param delta Convergence tolerance for the iterative algorithm (default = 1e-4)
#' @param maxiter Maximum number of iterations (default = 30)
#' @param C Numeric value for the constraint bound to ensure stability (default = 5)
#' @param verbose Logical indicating whether to print iteration progress (default = TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item FF_hat: Estimated factor matrix (n x q)
#'   \item BB_hat: Estimated loading matrix (p x q)
#'   \item alpha_hat: Estimated intercept vector (p x 1)
#'   \item obj_loglik: Log-likelihood value
#'   \item obj_pen: Penalized objective function value
#'   \item index: Indices of variables with zero loadings (selected out)
#'   \item df_est: Estimated degrees of freedom
#'   \item iter: Number of iterations performed
#' }
#'
#' @details
#' The function implements an alternating minimization algorithm that iteratively updates
#' the factor matrix F and loading matrix B until convergence. Missing values are imputed
#' using column means (rounded for count/binary data). The algorithm includes identifiability
#' constraints and various penalty functions for inducing sparsity in the loading matrix.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' n <- 100; p <- 50
#' x <- matrix(rnorm(n*p), n, p)
#'
#' # Fit sparse GFM with lasso penalty
#' result <- sparseGFM(x, type = "continuous", q = 2, penalty = "lasso", lambda = 0.5)
#' }
#'
#' @export
sparseGFM <- function(x, type = c("continuous", "count", "binary"),
                      q = 2, penalty = c("lasso", "SCAD", "MCP", "group", "agroup", "gSCAD", "agSCAD", "gMCP", "agMCP", "alasso", "glasso", "aglasso"),
                      lambda = 1, gam = 1, tau = NULL, mat_sd = 1,
                      delta = 1e-4, maxiter = 30, C = 5, verbose=TRUE) {

  # Check input parameters
  type <- match.arg(type)
  penalty <- match.arg(penalty)

  # Set default tau values based on penalty type
  if(is.null(tau)) {
    if(penalty %in% c("SCAD", "gSCAD", "agSCAD")) {
      tau = 3.7
    } else if(penalty %in% c("MCP", "gMCP", "agMCP")) {
      tau = 3
    }
  }

  # Get data dimensions
  n <- nrow(x)
  p <- ncol(x)
  xx=x
  # Initialize parameters based on data type
  # 检查是否有全空的行或列
  if(any(rowSums(is.na(xx)) == ncol(xx))) {
    stop("Error: There are rows with all missing values. Please remove or impute before processing.")
  }
  if(any(colSums(is.na(xx)) == nrow(xx))) {
    stop("Error: There are columns with all missing values. Please remove or impute before processing.")
  }

  if (type == "continuous") {
    # 对连续型数据处理缺失值 - 使用列均值填充
    if(any(is.na(xx))) {
      for(j in 1:ncol(xx)) {
        if(any(is.na(xx[, j]))) {
          xx[is.na(xx[, j]), j] <- mean(xx[, j], na.rm = TRUE)
        }
      }
    }

    # For continuous data - Gaussian family
    ########  continuous  ########
    F_int <- {tmp <- suppressWarnings(
      tryCatch(GFM::gfm(list(xx), types = "gaussian", q = q,
                        algorithm = "AM", verbose = FALSE),
               error = function(e) NULL));
    if (is.null(tmp) || any(!is.finite(tmp$hB)))
      suppressWarnings(
        GFM::gfm(list(xx), types = "gaussian", q = q, verbose = FALSE))
    else tmp}

    family_type <- gaussian()

    # Initialize model parameters
    FF_hat <- F_int$hH
    BB_hat <- F_int$hB
    alpha_hat <- F_int$hmu
    BB_int <- BB_hat
  } else if (type == "count") {
    # 对计数型数据处理缺失值 - 使用列均值四舍五入填充
    if(any(is.na(xx))) {
      for(j in 1:ncol(xx)) {
        if(any(is.na(xx[, j]))) {
          xx[is.na(xx[, j]), j] <- round(mean(xx[, j], na.rm = TRUE))
        }
      }
    }

    # For count data - Poisson family
    F_int <- {tmp <- suppressWarnings(
      tryCatch(GFM::gfm(list(xx), types = "poisson", q = q,
                        algorithm = "AM", verbose = FALSE),
               error = function(e) NULL));
    if (is.null(tmp) || any(!is.finite(tmp$hB)))
      suppressWarnings(
        GFM::gfm(list(xx), types = "poisson", q = q, verbose = FALSE))
    else tmp}

    family_type <- poisson()

    # Initialize model parameters
    FF_hat <- F_int$hH
    BB_hat <- F_int$hB
    alpha_hat <- F_int$hmu
    BB_int <- BB_hat
  } else if (type == "binary") {
    # 对二元数据处理缺失值 - 使用列均值四舍五入填充
    if(any(is.na(xx))) {
      for(j in 1:ncol(xx)) {
        if(any(is.na(xx[, j]))) {
          xx[is.na(xx[, j]), j] <- round(mean(xx[, j], na.rm = TRUE))
        }
      }
    }

    # For binary data - Binomial family

    ########  binary  ########
    F_int <- {tmp <- suppressWarnings(
      tryCatch(GFM::gfm(list(xx), types = "binomial", q = q,
                        algorithm = "AM", verbose = FALSE),
               error = function(e) NULL));
    if (is.null(tmp) || any(!is.finite(tmp$hB)))
      suppressWarnings(
        GFM::gfm(list(xx), types = "binomial", q = q, verbose = FALSE))
    else tmp}

    family_type <- binomial()

    # Initialize model parameters
    FF_hat <- F_int$hH
    BB_hat <- F_int$hB
    alpha_hat <- F_int$hmu
    BB_int <- BB_hat
  }

  # Ensure FF_hat is a matrix even if there's only one column
  if (!is.matrix(FF_hat)) {
    FF_hat <- as.matrix(FF_hat)
  }
  if (!is.matrix(BB_hat)) {
    BB_hat <- as.matrix(BB_hat)
  }
  #print(BB_hat)
  # Iterative optimization
  iter <- 0
  current_delta <- Inf
  while (current_delta > delta && iter < maxiter && !all(abs(BB_hat) < 1e-6)) {
    iter <- iter + 1
    #cat("iter = ", iter, "\n")
    # Store previous parameter values
    temp_FF <- FF_hat
    temp_BB <- BB_hat

    # Update B and alpha
    # for (j in 1:p) {
    #   grre <- glm(x[, j] ~ FF_hat, family = family_type)
    #
    #   # Check for NaN, Inf, NA in coefficients
    #   if(any(is.nan(coef(grre))) || any(is.infinite(coef(grre))) || any(is.na(coef(grre)))) {
    #     BB_hat[j, ] <- BB_hat[j, ] + 0.01
    #     alpha_hat[j] <- alpha_hat[j] + 0.01
    #   } else {
    #     # Projection and shrinkage processing
    #     coefs <- coef(grre)
    #     coef_norm <- sqrt(sum(coefs^2))
    #
    #     if(coef_norm > C) {
    #       coefs <- coefs * (C / coef_norm)
    #     }
    #
    #     BB_hat[j, ] <- coefs[-1]
    #     alpha_hat[j] <- coefs[1]
    #   }
    # }
    #
    #
    for (j in 1:p) {
      grre <- tryCatch(
        glm(x[, j] ~ FF_hat, family = family_type),
        error = function(e) NULL
      )

      if(is.null(grre) || any(is.nan(coef(grre))) || any(is.infinite(coef(grre))) || any(is.na(coef(grre)))) {
        BB_hat[j, ] <- BB_hat[j, ] + 0.01
        alpha_hat[j] <- alpha_hat[j] + 0.01
      } else {
        coefs <- coef(grre)
        coef_norm <- sqrt(sum(coefs^2))

        if(coef_norm > C) {
          coefs <- coefs * (C / coef_norm)
        }

        BB_hat[j, ] <- coefs[-1]
        alpha_hat[j] <- coefs[1]
      }
    }
    # cat("BB_hat")
    #  print(BB_hat)
    # Apply thresholding to B based on penalty
    if (penalty == "lasso") {
      # B <- matrix(0, nrow = nrow(BB_hat), ncol = ncol(BB_hat))
      # B[abs(BB_hat) > lambda] <- BB_hat[abs(BB_hat) > lambda] - sign(BB_hat[abs(BB_hat) > lambda]) * lambda
      # B <- BB_hat
      B=BB_hat
      B[abs(BB_hat) <= lambda] <- 0
      B[BB_hat > lambda] <- B[BB_hat > lambda] - lambda
      B[BB_hat < -lambda] <- B[BB_hat < -lambda] + lambda
      penalty_value <- lambda * sum(abs(B))
      BB_hat <- B
    } else if (penalty == "alasso") {
      ada_lambda=lambda/((abs(BB_int))^gam)
      B <- BB_hat
      B[abs(BB_hat) <= ada_lambda] <- 0
      B[BB_hat > ada_lambda] <- B[BB_hat > ada_lambda] - ada_lambda
      B[BB_hat < -ada_lambda] <- B[BB_hat < -ada_lambda] + ada_lambda
      penalty_value <- sum(ada_lambda * abs(B))
      BB_hat <- B
    } else if (penalty == "SCAD") {
      B <- BB_hat
      B[abs(BB_hat) <= lambda] <- 0
      B[BB_hat > lambda & BB_hat <= 2*lambda] <- B[BB_hat > lambda & BB_hat <= 2*lambda] - lambda
      B[BB_hat < -lambda & BB_hat >= -2*lambda] <- B[BB_hat < -lambda & BB_hat >= -2*lambda] + lambda

      idx <- abs(BB_hat) > 2*lambda & abs(BB_hat) <= tau*lambda
      if(any(idx)) {
        B[idx] <- sign(B[idx]) * ((tau-1) * abs(B[idx]) - tau*lambda) / (tau-2)
      }

      ab_vec <- abs(B)
      tem0_vec <- (lambda * ab_vec) * (ab_vec < lambda)
      tem1_vec <- ((tau*lambda*(ab_vec-lambda)-(ab_vec^2-lambda^2)/2)/(tau-1)+lambda^2) * (ab_vec >= lambda) * (ab_vec < tau*lambda)
      tem2_vec <- ((tau+1)*lambda^2/2) * (ab_vec >= tau*lambda)
      penalty_value <- sum(tem0_vec + tem1_vec + tem2_vec)

      BB_hat <- B
    } else if (penalty == "MCP") {
      B <- BB_hat
      B[abs(BB_hat) <= lambda] <- 0

      idx <- abs(BB_hat) > lambda & abs(BB_hat) <= tau*lambda
      if(any(idx)) {
        B[idx] <- sign(B[idx]) * (abs(B[idx]) - lambda) / (1 - 1/tau)
      }

      ab_vec <- abs(B)
      tem0_vec <- (lambda * ab_vec - ab_vec^2 / (2*tau)) * (ab_vec < tau*lambda)
      tem1_vec <- (tau * lambda^2 / 2) * (ab_vec >= tau*lambda)
      penalty_value <- sum(tem0_vec + tem1_vec)
      BB_hat <- B
    } else if (penalty %in% c("group", "glasso")) {
      row_12norm <- sqrt(rowSums(BB_hat^2))
      vec <- 1 - lambda / row_12norm
      vec[vec < 0] <- 0
      vec <- matrix(vec, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- BB_hat * vec
      penalty_value <- lambda * sum(sqrt(rowSums(B^2)))

      BB_hat <- B
    } else if (penalty %in% c("agroup", "aglasso")) {
      row_12norm <- sqrt(rowSums(BB_hat^2))
      wei <- sqrt(rowSums(as.matrix(BB_int)^2))^(-gam)
      vec <- 1 - ((wei*lambda) / row_12norm)
      vec[vec < 0] <- 0
      vec <- matrix(vec, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- BB_hat * vec
      penalty_value <- sum(lambda * wei * sqrt(rowSums(B^2)))

      BB_hat <- B
    } else if (penalty == "gSCAD") {
      B <- BB_hat
      row_12norm <- sqrt(rowSums(BB_hat^2))
      scad_row <- row_12norm

      scad_row[row_12norm <= lambda] <- 0
      scad_row[row_12norm > lambda & row_12norm <= 2*lambda] <- scad_row[row_12norm > lambda & row_12norm <= 2*lambda] - lambda

      idx <- row_12norm > 2*lambda & row_12norm <= tau*lambda
      if(any(idx)) {
        scad_row[idx] <- ((tau-1) * scad_row[idx] - tau*lambda) / (tau-2)
      }

      scal <- rep(1, length(row_12norm))
      scal[row_12norm > 0] <- scad_row[row_12norm > 0] / row_12norm[row_12norm > 0]
      scal_mat <- matrix(scal, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- scal_mat * BB_hat

      ab_vec <- abs(scad_row)
      tem0_vec <- (lambda * ab_vec) * (ab_vec < lambda)
      tem1_vec <- ((tau*lambda*(ab_vec-lambda)-(ab_vec^2-lambda^2)/2)/(tau-1)+lambda^2) * (ab_vec >= lambda) * (ab_vec < tau*lambda)
      tem2_vec <- ((tau+1)*lambda^2/2) * (ab_vec >= tau*lambda)
      penalty_value <- sum(tem0_vec + tem1_vec + tem2_vec)

      BB_hat <- B
    } else if (penalty == "agSCAD") {
      B <- BB_hat
      row_12norm <- sqrt(rowSums(BB_hat^2))
      wei <- sqrt(rowSums(as.matrix(BB_int)^2))^(-gam)
      ada_lambda <- wei * lambda
      scad_row <- row_12norm

      scad_row[row_12norm <= ada_lambda] <- 0
      idx1 <- row_12norm > ada_lambda & row_12norm <= 2*ada_lambda
      if(any(idx1)) {
        scad_row[idx1] <- scad_row[idx1] - ada_lambda[idx1]
      }

      idx2 <- row_12norm > 2*ada_lambda & row_12norm <= tau*ada_lambda
      if(any(idx2)) {
        scad_row[idx2] <- ((tau-1) * scad_row[idx2] - tau*ada_lambda[idx2]) / (tau-2)
      }

      scal <- rep(1, length(row_12norm))
      scal[row_12norm > 0] <- scad_row[row_12norm > 0] / row_12norm[row_12norm > 0]
      scal_mat <- matrix(scal, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- scal_mat * BB_hat

      ab_vec <- abs(scad_row)
      tem0_vec <- (ada_lambda * ab_vec) * (ab_vec < ada_lambda)
      tem1_vec <- ((tau*ada_lambda*(ab_vec-ada_lambda)-(ab_vec^2-ada_lambda^2)/2)/(tau-1)+ada_lambda^2) * (ab_vec >= ada_lambda) * (ab_vec < tau*ada_lambda)
      tem2_vec <- ((tau+1)*ada_lambda^2/2) * (ab_vec >= tau*ada_lambda)
      penalty_value <- sum(tem0_vec + tem1_vec + tem2_vec)

      BB_hat <- B
    } else if (penalty == "gMCP") {
      B <- BB_hat
      row_12norm <- sqrt(rowSums(BB_hat^2))
      mcp_row <- row_12norm

      mcp_row[row_12norm <= lambda] <- 0

      idx <- row_12norm > lambda & row_12norm <= tau*lambda
      if(any(idx)) {
        mcp_row[idx] <- (mcp_row[idx] - lambda) / (1-1/tau)
      }

      scal_mcp <- rep(1, length(row_12norm))
      scal_mcp[row_12norm > 0] <- mcp_row[row_12norm > 0] / row_12norm[row_12norm > 0]
      scal_mcp_mat <- matrix(scal_mcp, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- scal_mcp_mat * BB_hat

      ab_vec <- abs(mcp_row)
      tem0_vec <- (lambda * ab_vec - ab_vec^2 / (2*tau)) * (ab_vec < tau*lambda)
      tem1_vec <- (tau * lambda^2 / 2) * (ab_vec >= tau*lambda)
      penalty_value <- sum(tem0_vec + tem1_vec)

      BB_hat <- B
    } else if (penalty == "agMCP") {
      B <- BB_hat
      row_12norm <- sqrt(rowSums(BB_hat^2))
      wei <- sqrt(rowSums(as.matrix(BB_int)^2))^(-gam)
      ada_lambda <- wei * lambda
      mcp_row <- row_12norm

      mcp_row[row_12norm <= ada_lambda] <- 0

      idx <- row_12norm > ada_lambda & row_12norm <= tau*ada_lambda
      if(any(idx)) {
        mcp_row[idx] <- (mcp_row[idx] - ada_lambda[idx]) / (1-1/tau)
      }

      scal_mcp <- rep(1, length(row_12norm))
      scal_mcp[row_12norm > 0] <- mcp_row[row_12norm > 0] / row_12norm[row_12norm > 0]
      scal_mcp_mat <- matrix(scal_mcp, nrow = nrow(BB_hat), ncol = ncol(BB_hat), byrow = FALSE)
      B <- scal_mcp_mat * BB_hat

      ab_vec <- abs(mcp_row)
      tem0_vec <- (ada_lambda * ab_vec - ab_vec^2 / (2*tau)) * (ab_vec < tau*ada_lambda)
      tem1_vec <- (tau * ada_lambda^2 / 2) * (ab_vec >= tau*ada_lambda)
      penalty_value <- sum(tem0_vec + tem1_vec)

      BB_hat <- B
    }

    # Apply identifiability constraint if all elements of BB_hat are not very small
    # if (!all(abs(BB_hat) < 1e-6)) {
    #   id_constraint <- add_identifiability(FF_hat, BB_hat, alpha_hat)
    #   BB_hat <- id_constraint$B
    #   alpha_hat <- id_constraint$mu
    # }
    #cat("after id")
    #print(id_constraint)
    # Small values to zero
    #BB_hat[abs(BB_hat) < 1e-6] <- 0

    # Calculate delta for B
    delta_B <- base::norm((BB_hat - temp_BB), type = "F") / max(base::norm(temp_BB, type = "F"), 1e-10)

    # Update F
    # for (i in 1:n) {
    #   tryCatch({
    #     mod_f <- glm.fit(BB_hat, x[i, ], family = family_type, intercept = FALSE, offset = alpha_hat)
    #     coefs <- coef(mod_f)
    #
    #     # Check for NaN, Inf, NA in coefficients
    #     if(any(is.nan(coefs)) || any(is.infinite(coefs)) || any(is.na(coefs))) {
    #       FF_hat[i, ] <- FF_hat[i, ] + 0.01
    #     } else {
    #       # Projection and shrinkage processing
    #       coef_norm <- sqrt(sum(coefs^2)+1)
    #
    #       if(coef_norm > C) {
    #         coefs <- coefs * (C / coef_norm)
    #       }
    #
    #       FF_hat[i, ] <- coefs
    #     }
    #   }, error = function(e) {
    #     FF_hat[i, ] <- FF_hat[i, ] + 0.01
    #   })
    # }
    for (i in 1:n) {
      mod_f <- tryCatch(
        glm.fit(BB_hat, x[i, ], family = family_type, intercept = FALSE, offset = alpha_hat),
        error = function(e) NULL
      )

      if(is.null(mod_f) || any(is.nan(coef(mod_f))) || any(is.infinite(coef(mod_f))) || any(is.na(coef(mod_f)))) {
        FF_hat[i, ] <- FF_hat[i, ] + 0.01
      } else {
        coefs <- coef(mod_f)
        coef_norm <- sqrt(sum(coefs^2)+1)

        if(coef_norm > C) {
          coefs <- coefs * (C / coef_norm)
        }

        FF_hat[i, ] <- coefs
      }
    }

    # Apply identifiability constraint to F if all elements of BB_hat are not very small
    if (!all(abs(BB_hat) < 1e-6)) {
      FF_hat <- add_identifiability(FF_hat, BB_hat, alpha_hat)$H
    }
    #cat("after id")
    #print(id_constraint$H)
    # Calculate delta for F
    delta_F <- base::norm((FF_hat - temp_FF), type = "F") / max(base::norm(temp_FF, type = "F"), 1e-10)

    # Calculate overall delta
    current_delta <- delta_B + delta_F
    if(verbose == TRUE){
      cat("Current delta:", current_delta, "\n")
    }
  }

  # Calculate parameter estimates
  mat_para_est <- FF_hat %*% t(BB_hat) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha_hat))

  # Calculate log-likelihood
  if (type == "continuous") {
    obj_loglik <- sum(log(dnorm(x, mat_para_est, sd = mat_sd) + 1e-8),na.rm = T)
  } else if (type == "count") {
    obj_loglik <- sum(log(dpois(x, exp(mat_para_est)) + 1e-8),na.rm = T)
  } else if (type == "binary") {
    obj_loglik <- sum(log(dbinom(x, 1, plogis(mat_para_est)) + 1e-8),na.rm = T)
  }

  # Calculate objective function
  obj_pen <- obj_loglik - n * penalty_value

  # Get indices of zero rows in BB_hat
  row_norms <- sqrt(rowSums(BB_hat^2))
  index <- which(row_norms < 1e-6)

  df_est=sum(sqrt(rowSums(BB_hat^2))!= 0)+sum((sqrt(rowSums(BB_hat^2))/sqrt(rowSums(BB_int^2)))*(q-1)) + p+n*q


  # Return results
  return(list(
    FF_hat = FF_hat,
    BB_hat = BB_hat,
    alpha_hat = alpha_hat,
    obj_loglik = obj_loglik,
    obj_pen = obj_pen,
    index = index,
    df_est = df_est,
    iter = iter
  ))
}


#' Cross-Validation for Sparse Generalized Factor Model
#'
#' @description
#' Performs cross-validation to select the optimal lambda parameter for sparse generalized
#' factor models using SIC (Sparsity Information Criterion).
#'
#' @param x A numeric matrix of observations (n x p)
#' @param type Character string specifying the data type ("continuous", "count", or "binary")
#' @param q Integer specifying the number of latent factors (default = 2)
#' @param penalty Character string specifying the penalty type (see sparseGFM for options)
#' @param lambda_range Numeric vector of lambda values to evaluate (default = seq(0.1, 1, by = 0.1))
#' @param gam Numeric value for the adaptive weight parameter (default = 1)
#' @param tau Numeric value for the shape parameter in SCAD/MCP penalties
#' @param mat_sd Standard deviation for continuous data (default = 1)
#' @param delta Convergence tolerance (default = 1e-4)
#' @param maxiter Maximum number of iterations (default = 30)
#' @param C Constraint bound for stability (default = 5)
#' @param bic_type Character string specifying BIC type. Options are:
#'        \itemize{
#'          \item "dd": Default SIC using estimated degrees of freedom
#'          \item "as": Alternative simplified BIC
#'        }
#' @param verbose Logical indicating whether to print progress (default = TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item optimal_lambda: Selected optimal lambda value
#'   \item optimal_model: Model fitted with optimal lambda
#'   \item all_results: List of all fitted models for each lambda
#'   \item objloglik: Vector of log-likelihood values for each lambda
#'   \item bic_dd: Vector of default SIC values
#'   \item bic_as: Vector of alternative SIC values
#'   \item df_dd: Vector of degrees of freedom (default method)
#'   \item df_as: Vector of degrees of freedom (alternative method)
#'   \item lambda_range: Vector of evaluated lambda values
#' }
#'
#' @details
#' The function fits sparse GFM models for each lambda value in lambda_range and
#' calculates two types of SIC for model selection. The optimal lambda is chosen
#' as the one minimizing the selected SIC criterion.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' n <- 100; p <- 50
#' x <- matrix(rnorm(n*p), n, p)
#'
#' # Cross-validation for lambda selection
#' cv_result <- cv.sparseGFM(x, type = "continuous", q = 2, penalty = "lasso")
#' optimal_lambda <- cv_result$optimal_lambda
#' }
#'
#' @export
cv.sparseGFM <- function(x,
                         type = c("continuous", "count", "binary"),
                         q = 2,
                         penalty = c("lasso", "SCAD", "MCP", "group", "agroup", "gSCAD", "agSCAD", "gMCP", "agMCP", "alasso", "glasso", "aglasso"),
                         lambda_range = seq(0.1, 1, by = 0.1),
                         gam = 1,
                         tau = NULL,
                         mat_sd = 1,
                         delta = 1e-4,
                         maxiter = 30,
                         C = 5,
                         bic_type = "dd", verbose=TRUE) {

  # Match type and penalty arguments
  type <- match.arg(type)
  penalty <- match.arg(penalty)

  # Get dimensions
  n <- nrow(x)
  p <- ncol(x)

  # Initialize vectors to store BIC values
  bic_dd <- numeric(length(lambda_range))
  bic_as <- numeric(length(lambda_range))
  df_dd <- numeric(length(lambda_range))
  df_as <- numeric(length(lambda_range))
  objloglik=numeric(length(lambda_range))
  # Initialize list to store all estimation results
  all_results <- list()

  # Loop through each lambda value
  for (i in seq_along(lambda_range)) {
    # Get current lambda
    lambda <- lambda_range[i]

    # Call sparseGFM with current lambda
    result <- sparseGFM(x = x,
                        type = type,
                        q = q,
                        penalty = penalty,
                        lambda = lambda,
                        gam = gam,
                        tau = tau,
                        mat_sd = mat_sd,
                        delta = delta,
                        maxiter = maxiter,
                        C = C,verbose = verbose)

    # Store the result
    all_results[[i]] <- result

    # Calculate BIC values
    # Default BIC (bic_dd)
    objloglik[i] = result$obj_loglik/(n*p)
    s_hat=(p-length(result$index))
    bic_dd[i] <- -2 * result$obj_loglik/(n*p) + result$df_est * log(n*s_hat/(n+p))/(n*p)
    df_dd[i] = result$df_est
    # Alternative BIC (bic_as)
    #cat("length(result$index) = ",length(result$index), "\n")

    df_as[i] =  (p - length(result$index) + n*q + p*q)

    #cat("df_as[i]",df_as[i] ,"\n")
    bic_as[i] <- -2 * result$obj_loglik/(n*p) +  (p - length(result$index) + n*q + p*q) * log(n*s_hat/(n+p))/(n*p)
    # Print progress
    if (verbose==TRUE){
      cat("Lambda:", lambda, "- BIC_dd:", bic_dd[i],
          "- BIC_as:", bic_as[i], "\n")
    }
  }

  # Determine which BIC to use for selection
  if (bic_type == "as") {
    optimal_index <- which.min(bic_as)
    optimal_bic <- bic_as[optimal_index]
    all_bic <- bic_as
  } else {
    optimal_index <- which.min(bic_dd)
    optimal_bic <- bic_dd[optimal_index]
    all_bic <- bic_dd
  }

  optimal_lambda <- lambda_range[optimal_index]
  optimal_model <- all_results[[optimal_index]]

  # Return results
  return(list(
    optimal_lambda = optimal_lambda,
    optimal_model = optimal_model,
    all_results = all_results,
    objloglik = objloglik,
    bic_dd = bic_dd,
    bic_as = bic_as,
    df_dd = df_dd,
    df_as = df_as,
    lambda_range = lambda_range
  ))
}



#' Factor Number Selection for Sparse Generalized Factor Model
#'
#' @description
#' Determines the optimal number of factors for sparse generalized factor models
#' using information criteria (SIC - Sparsity Information Criterion variants).
#'
#' @param x A numeric matrix of observations (n x p)
#' @param type Character string specifying the data type ("continuous", "count", or "binary")
#' @param q_range Integer vector of factor numbers to evaluate (default = 1:5)
#' @param penalty Character string specifying the penalty type (see sparseGFM for options)
#' @param lambda_range Numeric vector of lambda values for cross-validation (default = seq(0.1, 1, by = 0.1))
#' @param gam Numeric value for the adaptive weight parameter (default = 1)
#' @param tau Numeric value for the shape parameter in SCAD/MCP penalties
#' @param mat_sd Standard deviation for continuous data (default = 1)
#' @param delta Convergence tolerance (default = 1e-4)
#' @param maxiter Maximum number of iterations (default = 30)
#' @param C Constraint bound for stability (default = 5)
#' @param sic_type Character string specifying SIC type for selection. Options are:
#'        \itemize{
#'          \item "sic1": SIC with estimated df and (n+p) denominator
#'          \item "sic2": SIC with simplified df and (n+p) denominator
#'          \item "sic3": SIC with estimated df and max(n,p) denominator
#'          \item "sic4": SIC with simplified df and max(n,p) denominator
#'        }
#' @param verbose Logical indicating whether to print progress (default = TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item optimal_q: Selected optimal number of factors
#'   \item optimal_model: Model fitted with optimal q
#'   \item all_results: List of all fitted models for each q
#'   \item objpen: Vector of penalized objective values
#'   \item objlogli: Vector of log-likelihood values
#'   \item sic1, sic2, sic3, sic4: Vectors of different SIC values
#'   \item sic21, sic22, sic23, sic24: Alternative SIC formulations
#'   \item lambda_opt: Vector of optimal lambda values for each q
#'   \item df_dd: Vector of degrees of freedom (default method)
#'   \item df_as: Vector of degrees of freedom (alternative method)
#'   \item q_range: Vector of evaluated q values
#'   \item used_sic_type: The SIC type used for selection
#'   \item optimal_sic_value: The optimal SIC value
#' }
#'
#' @details
#' For each q value, the function performs cross-validation to select optimal lambda,
#' then calculates various SIC measures. The optimal q minimizes the selected SIC.
#' This provides automatic selection of the latent dimension in factor models.
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' n <- 100; p <- 50
#' x <- matrix(rnorm(n*p), n, p)
#'
#' # Select optimal number of factors
#' facnum_result <- facnum.sparseGFM(x, type = "continuous", q_range = 1:5, penalty = "lasso")
#' optimal_q <- facnum_result$optimal_q
#' }
#'
#' @export
facnum.sparseGFM <- function(x,
                             type = c("continuous", "count", "binary"),
                             q_range = 1:5,
                             penalty = c("lasso", "SCAD", "MCP", "group", "agroup", "gSCAD", "agSCAD", "gMCP", "agMCP", "alasso", "glasso", "aglasso"),
                             lambda_range = seq(0.1, 1, by = 0.1),
                             gam = 1,
                             tau = NULL,
                             mat_sd = 1,
                             delta = 1e-4,
                             maxiter = 30,
                             C = 5,
                             sic_type = "sic1",
                             verbose=TRUE) {

  # Match type and penalty arguments
  type <- match.arg(type)
  penalty <- match.arg(penalty)

  # Get dimensions
  n <- nrow(x)
  p <- ncol(x)

  # Initialize vectors to store SIC values
  sic1 <- numeric(length(q_range))
  sic2 <- numeric(length(q_range))
  sic3 <- numeric(length(q_range))
  sic4 <- numeric(length(q_range))
  sic21 <- numeric(length(q_range))
  sic22 <- numeric(length(q_range))
  sic23 <- numeric(length(q_range))
  sic24 <- numeric(length(q_range))
  objpen = numeric(length(q_range))
  objlogli = numeric(length(q_range))
  df_dd = numeric(length(q_range))
  df_as = numeric(length(q_range))
  lambda_opt = numeric(length(q_range))
  # Initialize list to store all estimation results
  all_results <- list()

  # Loop through each q value
  for (i in seq_along(q_range)) {
    # Get current q
    q <- q_range[i]

    # Call sparseGFM with current q
    # result <- sparseGFM(x = x,
    #                     type = type,
    #                     q = q,
    #                     penalty = penalty,
    #                     lambda = lambda,
    #                     gam = gam,
    #                     tau = tau,
    #                     mat_sd = mat_sd,
    #                     delta = delta,
    #                     maxiter = maxiter,
    #                     C = C,
    #                     verbose = verbose)
    cv_result <- cv.sparseGFM(x,
                              type = type,
                              q = q,
                              penalty = penalty,
                              lambda_range = lambda_range,
                              gam = gam,
                              tau = tau,
                              mat_sd = mat_sd,
                              delta = delta,
                              maxiter = maxiter,
                              C = C,
                              verbose = verbose)
    result = cv_result$optimal_model
    lambda_opt[i]=cv_result$optimal_lambda
    index_model = which(lambda_range == lambda_opt[i])
    df_dd[i] = (cv_result$df_dd)[index_model]
    df_as[i] = (cv_result$df_as)[index_model]
    #cat("df_dd" , df_dd[i], "\n")
    #cat("df_as" , df_as[i], "\n")
    # Store the result
    all_results[[i]] <- result

    # Calculate SIC values
    objpen[i]=result$obj_pen/(n*p)
    objlogli[i]=result$obj_loglik/(n*p)
    s_hat=(p-length(result$index))
    # sic1[i] <- -2 * result$obj_pen/(n*p) + result$df_est * log(n*p/(n+p))/(n*p)
    # sic2[i] <- -2 * result$obj_pen/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*p/(n+p))/(n*p)
    # sic3[i] <- -2 * result$obj_pen/(n*p) + result$df_est * log(n*p/(max(n,p)))/(n*p)
    # sic4[i] <- -2 * result$obj_pen/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*p/(max(n,p)))/(n*p)
    #
    sic1[i] <- -2 * result$obj_loglik/(n*p) + result$df_est * log(n*s_hat/(n+p))/(n*p)
    sic2[i] <- -2 * result$obj_loglik/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*s_hat/(n+p))/(n*p)
    sic3[i] <- -2 * result$obj_loglik/(n*p) + result$df_est * log(n*s_hat/(max(n,p)))/(n*p)
    sic4[i] <- -2 * result$obj_loglik/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*s_hat/(max(n,p)))/(n*p)



    sic21[i] <- -result$obj_loglik/(n*p) + result$df_est * log(n*s_hat/(n+p))/(n*p)
    sic22[i] <- -result$obj_loglik/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*s_hat/(n+p))/(n*p)
    sic23[i] <- -result$obj_loglik/(n*p) + result$df_est * log(n*s_hat/(max(n,p)))/(n*p)
    sic24[i] <- -result$obj_loglik/(n*p) + (p - length(result$index) + n*q + p*q) * log(n*s_hat/(max(n,p)))/(n*p)

    if (verbose==TRUE){
      # Print progress
      cat("q =", q, "\n")
      cat("SIC1:", sic1[i], "\n")
      cat("SIC2:", sic2[i], "\n")
      cat("SIC3:", sic3[i], "\n")
      cat("SIC4:", sic4[i], "\n")
      cat("-----------------------\n")
    }
  }

  # Determine which SIC to use for selection
  if (sic_type == "sic1") {
    optimal_index <- which.min(sic1)
    optimal_sic <- sic1[optimal_index]
    used_sic <- sic1
  } else if (sic_type == "sic2") {
    optimal_index <- which.min(sic2)
    optimal_sic <- sic2[optimal_index]
    used_sic <- sic2
  } else if (sic_type == "sic3") {
    optimal_index <- which.min(sic3)
    optimal_sic <- sic3[optimal_index]
    used_sic <- sic3
  } else if (sic_type == "sic4") {
    optimal_index <- which.min(sic4)
    optimal_sic <- sic4[optimal_index]
    used_sic <- sic4
  }

  optimal_q <- q_range[optimal_index]
  optimal_model <- all_results[[optimal_index]]

  # Return results
  return(list(
    optimal_q = optimal_q,
    optimal_model = optimal_model,
    all_results = all_results,
    objpen = objpen,
    objlogli = objlogli,
    sic1 = sic1,
    sic2 = sic2,
    sic3 = sic3,
    sic4 = sic4,
    sic21 = sic21,
    sic22 = sic22,
    sic23 = sic23,
    sic24 = sic24,
    lambda_opt = lambda_opt,
    df_dd = df_dd,
    df_as = df_as,
    q_range = q_range,
    used_sic_type = sic_type,
    optimal_sic_value = optimal_sic
  ))
}







#' Apply Identifiability Constraints
#' @description Internal function to ensure identifiability of factor and loading matrices
#' @param H Factor matrix (n x q)
#' @param B Loading matrix (p x q)
#' @param mu Intercept vector (p x 1)
#' @return List with constrained H, B, and mu matrices


add_identifiability <- function(H, B, mu){
  # Load the irlba library
  #library(irlba)

  # Perform SVD decomposition with rank k = 10

  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  signB1 <- sign(svdHB$v[1,])
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)

  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)

  return(list(H=H, B=B, mu=mu))
}

#' Create Diagonal Matrix
#' @description Internal helper to create diagonal matrix from vector
#' @param vec Numeric vector
#' @return Diagonal matrix with vec on diagonal
#' @noRd
#' @keywords internal
Diag=function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}

#' Measure Factor Space Recovery
#' @description Calculate recovery metrics between estimated and true factor spaces
#' @param hH Estimated factor matrix
#' @param H True factor matrix
#' @param type Metric type ("trace_statistic" or "ccor")
#' @return Numeric recovery measure
#' @noRd
#' @keywords internal
measurefun <- function(hH, H, type=c('trace_statistic','ccor')){
  type <- match.arg(type)
  q <- min(ncol(H), ncol(hH))
  switch(type,
         ccor=cancor(hH, H)$cor[q],
         trace_statistic= trace_statistic_fun(hH, H))
}

#' Trace Statistic Function
#' @description Calculate trace statistic for factor space comparison
#' @param H Estimated factor matrix
#' @param H0 True factor matrix
#' @return Trace statistic value
#' @noRd
#' @keywords internal
trace_statistic_fun <- function(H, H0){

  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% ginv(t(H) %*% H) %*% t(H) %*% H0

  tr_fun(mat1) / tr_fun(t(H0) %*% H0)

}

#' Evaluate Subspace Distance
#' @description Calculate multiple distance metrics between two subspaces
#' @param A First subspace matrix
#' @param B Second subspace matrix
#' @param orthnm Whether to orthonormalize matrices (default TRUE)
#' @return Vector of distance metrics
#' @noRd
#' @keywords internal
eval.space <- function(A, B, orthnm = TRUE)
{
  if(!is.matrix(A)) A <- as.matrix(A)
  if(!is.matrix(B)) B <- as.matrix(B)

  if(orthnm)
  {
    A <- qr.Q(qr(A))
    B <- qr.Q(qr(B))
  }

  mat <- t(B) %*% A %*% t(A) %*% B
  d <- eigen(mat)$values
  d <- (d + abs(d))/2
  q <- sqrt(prod(d))
  r <- sqrt(mean(d))

  c(q, r, acos(q))
}



#' Evaluate Binary Classification Performance
#' @description Calculate various performance metrics for variable selection
#' @param true_vector True selection indicator vector
#' @param pred_vector Predicted selection indicator vector
#' @return List containing confusion matrix and performance metrics (accuracy, precision, recall, etc.)

evaluate_performance <- function(true_vector, pred_vector) {
  # 确保向量是二进制值（0或1）
  true_binary <- as.integer(true_vector > 0)
  pred_binary <- as.integer(pred_vector > 0)

  # 计算混淆矩阵元素
  TP <- sum(true_binary == 1 & pred_binary == 1)  # 真正例
  TN <- sum(true_binary == 0 & pred_binary == 0)  # 真负例
  FP <- sum(true_binary == 0 & pred_binary == 1)  # 假正例
  FN <- sum(true_binary == 1 & pred_binary == 0)  # 假负例

  # 计算评估指标
  accuracy <- (TP + TN) / (TP + TN + FP + FN)  # 准确率

  # 防止除零错误
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), NA)  # 精确率
  recall_tpr <- ifelse(TP + FN > 0, TP / (TP + FN), NA)  # 召回率/真正例率(TPR)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), NA)  # 特异度
  fpr <- ifelse(FP + TN > 0, FP / (FP + TN), NA)  # 假正例率(FPR)

  # F1分数
  f1_score <- ifelse(precision + recall_tpr > 0,
                     2 * (precision * recall_tpr) / (precision + recall_tpr), NA)

  # MCC (Matthews Correlation Coefficient)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denominator > 0,
                (TP * TN - FP * FN) / denominator, NA)

  # 返回结果
  return(list(
    confusion_matrix = matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE,
                              dimnames = list(c("Actual Positive", "Actual Negative"),
                                              c("Predicted Positive", "Predicted Negative"))),
    accuracy = accuracy,
    precision = precision,
    recall_tpr = recall_tpr,
    specificity = specificity,
    fpr = fpr,
    f1_score = f1_score,
    mcc = mcc
  ))
}







