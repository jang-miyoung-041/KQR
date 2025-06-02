l2norm = function(x){
  return(sum(x^2))
}

gaussian_kernel = function(x1, x2, rho){
  return(exp(-rho * l2norm(x1 - x2)))
}

gaussian_kernel_mat = function(x, rho){
  if (class(x)[1] == "numeric")
    x = matrix(x)

  n_row = nrow(x)
  
  K = matrix(0, nrow = n_row, ncol = n_row)
  
  for (i in 1:n_row){
    for (j in 1:n_row){
      K[i,j] = gaussian_kernel(x[i,],x[j,], rho)
    }
  }
  return(K)
}

gaussian_kernel_mat2 = function(z, x, rho){
  if (class(x)[1] == "numeric")
    x = matrix(x)
  if (class(z)[1] == "numeric")
    z = matrix(z)
  
  
  K = matrix(0, nrow = nrow(z), ncol = nrow(x))
  
  for (i in 1:nrow(z)){
    for (j in 1:nrow(x)){
      K[i, j] = gaussian_kernel(z[i, ], x[j, ], rho)
    }
  }
  return(K)
}

#' @examples
#' set.seed(923)
#' #' n = 30
#' x_values = sort(runif(n, 0, 1))
#' y_values = sin(2 * pi * x_values) + cos(4 * pi * x_values) + rnorm(n, sd = 0.2)
#' # kernel specification
#' lambda = 0.001
#' rho = 1
#' # model fitting
#' model = fit_krr(x_values, y_values, lambda, rho)
#' print(model)
#'
#' @export
fit_krr = function(x, y, lambda, rho){
  
  x = matrix(x)
  
  n_row = nrow(x)
  
  K = gaussian_kernel_mat(x, rho)

  if (length(lambda) == 1){
    
    hat_alpha = solve(K + lambda * diag(n_row)) %*% y
    
    f_hat = t(hat_alpha) %*% K
    
    return(list("x" = x,
                "alpha_hat" = hat_alpha, 
                "best lambda" = lambda))
  } else{
    AIC = c()
    rss = c()
    hat_alpha = c()
    
    for (i in 1:length(lambda)){
      
      hat_alpha = solve(K + lambda[i] * diag(n_row)) %*% y
      
      f_hat = t(hat_alpha) %*% K
      
      rss[i] = sum((y - f_hat)^2)
      
      df = sum(diag(K %*% solve(K + lambda[i] * diag(n_row))))
      
      AIC[i] = n_row * log(rss[i] / n_row) + 2 * df
    }
    
    
    best_index = which.min(AIC)
    best_lambda = lambda[best_index]
    best_alpha = solve(K + best_lambda * diag(n_row)) %*% y
    
    
    return(list("x" = x,
                "alpha_hat" = best_alpha,
                "best_lambda" = best_lambda))
  }
}

#' new_x = seq(0, 1, len = 200)
#' prdict_k = predict_kernel(model, new_x, rho)
predict_kernel = function(model, new_x, rho){
  K.pred = gaussian_kernel_mat2(new_x, model$x, rho)
  f.pred = K.pred %*% model$alpha_hat
  
  return(f.pred)
}
 
#' grid_x = new_x
#' plot_spline(x_values, y_values, model, grid_x, rho)
plot_spline = function(x_values, y_values, model, grid_x, rho)
{
  y_pred = predict_kernel(model, grid_x, rho)
  data_plot = data.frame(x = x_values, y = y_values)
  spline_plot = data.frame(x = grid_x, y = y_pred)
  
  ggplot() +
    geom_point(data = data_plot, aes(x, y), color = "black") +
    geom_line(data = spline_plot, aes(x, y), color = "blue", linewidth = 1.2) +
    labs(title = "Fitted Kernel Regression", x = "x", y = "y") +
    xlim(c(min(x_values), max(x_values))) +
    theme_minimal()
}
