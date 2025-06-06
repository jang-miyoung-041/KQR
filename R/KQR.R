l2norm = function(x){
  return(sum(x^2))
}

gaussian_kernel = function(x1, x2, rho = 1){
  return(exp(-rho * l2norm(x1 - x2)))
}

poly_kernel = function(x1, x2, degree = 3){
  return((1 + x1 %*% t(x2))^degree)
}

rad_kernel = function(x1, x2, sigma2 = 1){
  return(exp(-l2norm(x1 - x2)/(2*sigma2)))
}

#' Generate the kernel matrix from the given data points using a specified kernel function
#'
#' This function computes the kernel matrix \(K\) where \(K_{ij} = \mathcal{K}(x_i, x_j)\) 
#' for a given set of input vectors using the specified kernel type.
#'
#' @param x A numeric vector or matrix containing the input data points. Each row is a data point.
#' @param kernel A character string specifying the kernel type. Options: "gaussian", "poly", or "rad".
#' @param rho A numeric value for the Gaussian kernel parameter (used when kernel = "gaussian").
#' @param degree An integer for the degree of the polynomial kernel (used when kernel = "poly").
#' @param sigma2 A numeric value for the variance in the radial basis kernel (used when kernel = "rad").
#'
#' @return A square symmetric kernel matrix of dimension n × n, where n is the number of data points.
#'
#' @examples
#' x = runif(50, 0, 1)
#' K = create_kernel_mat(x, kernel = "gaussian", rho = 0.5)
#' 
#' # For multivariate data:
#' x2 = matrix(runif(100), ncol = 2)
#' K2 = create_kernel_mat(x2, kernel = "poly", degree = 2)
#'
#' @export

create_kernel_mat = function(x, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1){
  if (class(x)[1] == "numeric")
    x = matrix(x)
  
  if (kernel == "gaussian") g = gaussian_kernel
  else if (kernel == "poly") g = poly_kernel
  else g = rad_kernel
  
  n_row = nrow(x)
  
  K = matrix(0, nrow = n_row, ncol = n_row)
  
  for (i in 1:n_row){
    for (j in 1:n_row){
      K[i,j] = g(x[i, ], x[j, ])
    }
  }
  return(K = K)
}

#' Check loss function for quantile regression
#'
#' @param y Actual values
#' @param y_hat Predicted values
#' @param tau Quantile level (0 < tau < 1)
#'
#' @return Numeric vector of check loss values
#'
#' @export

rho_tau = function(y, y_hat, tau = 0.6) {
  ifelse(y > y_hat, tau * abs(y - y_hat), -(1 - tau) * abs(y - y_hat))
}

#' Fit Kernel Quantile Regression via Linear Programming (No Intercept)
#'
#' @param K Kernel matrix (n x n)
#' @param y Response vector
#' @param tau Quantile level (e.g., 0.5)
#' @param lambda Regularization parameter (positive scalar)
#'
#' @return A list containing theta (dual coefficients), status, and objective value
#' 
#' @export
fit_kqr_lp = function(K, y, tau = 0.5, lambda = 1) {
  n <- length(y)
  library(lpSolve)
  ## Total variable is three
  total_var = 3 * n
  ## 0 (theta), tau (xi⁺), 1 - tau (xi⁻)
  obj = c(rep(0, n), rep(tau, n), rep(1 - tau, n))
  ## minimize function
  constr_mat = cbind(K / lambda, diag(n), -diag(n))  
  constr_dir = rep("=", n)
  constr_rhs = y
  # subject to
  lower_bounds = c(rep(-(1 - tau), n), rep(0, 2 * n))
  upper_bounds = c(rep(tau, n), rep(Inf, 2 * n))
  # LP
  result = lp("min",
              objective.in = obj,
              const.mat = constr_mat,
              const.dir = constr_dir,
              const.rhs = constr_rhs,
              transpose.constraints = FALSE)
  
  if (result$status != 0) warning("LP did not converge: status ", result$status)
  
  sol = result$solution
  list(
    theta = sol[1:n],
    objval = result$objval,
    status = result$status
  )
}

create_kernel_mat2 = function(z, x, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1){
  if (class(x)[1] == "numeric")
    x = matrix(x)
  if (class(z)[1] == "numeric")
    z = matrix(z)
  
  if (kernel == "gaussian") g = gaussian_kernel
  else if (kernel == "poly") g = poly_kernel
  else g = rad_kernel
  
  K = matrix(0, nrow = nrow(z), ncol = nrow(x))
  
  for (i in 1:nrow(z)){
    for (j in 1:nrow(x)){
      K[i, j] = gaussian_kernel(z[i, ], x[j, ], rho)
    }
  }
  return(K = K)
}

#' Predict using fitted Kernel Quantile Regression model
#'
#' @param x_new New data (matrix or numeric vector)
#' @param x_train Training data used to fit the model
#' @param theta Dual coefficients from KQR fit (fit$theta)
#' @param kernel Kernel type: "gaussian", "poly", or "rad"
#' @param rho Parameter for Gaussian kernel
#' @param degree Degree for polynomial kernel
#' @param sigma2 Variance for radial kernel
#' @param lambda Regularization parameter used in fitting
#'
#' @return A numeric vector of predicted values
#' 
#' @export
predict_kqr <- function(x_new, x_train, theta, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1, lambda = 1) {
  
  K_pred <- create_kernel_mat2(z = x_new, x = x_train,
                               kernel = kernel, rho = rho,
                               degree = degree, sigma2 = sigma2)
  
  f_pred <- (1 / lambda) * K_pred %*% theta
  
  return(as.vector(f_pred))
}

#' Plot Kernel Quantile Regression curves for multiple tau
#'
#' @param x_train Matrix of training inputs (n x p)
#' @param y_train Vector of responses (length n)
#' @param x_grid Points to predict (m x p)
#' @param taus Numeric vector of quantile levels (e.g., c(0.1, 0.5, 0.9))
#' @param kernel Kernel type ("gaussian", "poly", "rad")
#' @param rho Gaussian kernel parameter
#' @param degree Polynomial kernel degree
#' @param sigma2 Radial kernel variance
#' @param lambda Regularization parameter
#'
#' @return A ggplot object
#' 
#' @export
plot_kqr <- function(x_train, y_train, x_grid, taus = c(0.1, 0.5, 0.9),
                     kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1, lambda = 1) {
  plot_data <- data.frame()
  
  for (tau in taus) {
    K <- create_kernel_mat(x_train, kernel = kernel, rho = rho, degree = degree, sigma2 = sigma2)
    fit <- fit_kqr_no_intercept(K, y_train, tau = tau, lambda = lambda)
    
    K_pred <- create_kernel_mat2(x_grid, x_train, kernel = kernel, rho = rho, degree = degree, sigma2 = sigma2)
    y_hat <- (1 / lambda) * K_pred %*% fit$theta
    
    plot_data <- rbind(plot_data,
                       data.frame(x = x_grid[, 1], y = as.vector(y_hat), tau = factor(tau)))
  }
  
  df_raw <- data.frame(x = x_train[, 1], y = y_train)
  
  ggplot() +
    geom_point(data = df_raw, aes(x = x, y = y), color = "black", alpha = 0.6) +
    geom_line(data = plot_data, aes(x = x, y = y, color = tau), size = 1) +
    labs(title = "Kernel Quantile Regression Curves", x = "x", y = "f(x)", color = "τ") +
    theme_minimal()
}


