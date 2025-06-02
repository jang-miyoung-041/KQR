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

create_kernel_mat = function(x, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1){
  if (class(x)[1] == "numeric")
    x = matrix(x)
  
  if (kernel == "gaussian") g = gaussian_kernel
  else if (kernel == "ploy") g = poly_kernel
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

rho_tau = function(y, y_hat, tau = 0.6) {
  ifelse(y > y_hat, tau * abs(y - y_hat), -(1 - tau) * abs(y - y_hat))
}

predict_kqr <- function(alpha, x_train, x_new, kernel = "gaussian", param = list(rho = 1)) {
  n_train <- nrow(x_train)
  n_new <- nrow(x_new)
  f_pred <- numeric(n_new)
  
  for (i in 1:n_new) {
    total <- 0
    for (j in 1:n_train) {
      k_val <- switch(kernel,
                      "gaussian" = exp(-param$rho * sum((x_new[i, ] - x_train[j, ])^2)),
                      "poly"     = (1 + sum(x_new[i, ] * x_train[j, ]))^param$degree,
                      "rbf"      = exp(-sum((x_new[i, ] - x_train[j, ])^2) / (2 * param$sigma2)),
                      stop("Unknown kernel")
      )
      total <- total + alpha[j] * k_val
    }
    f_pred[i] <- total
  }
  
  return(f_pred)
}