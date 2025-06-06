\name{predict_kqr}
\alias{predict_kqr}
\title{Predict using fitted Kernel Quantile Regression model}
\usage{
predict_kqr(x_new, x_train, theta, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1, lambda = 1)
}
\arguments{
  \item{x_new}{A numeric vector or matrix of new input data for prediction. Each row is an observation.}
  \item{x_train}{A numeric matrix or vector used in training to compute the kernel matrix.}
  \item{theta}{A numeric vector of dual coefficients estimated from KQR fitting.}
  \item{kernel}{Kernel type to be used: one of \code{"gaussian"}, \code{"poly"}, or \code{"rad"}.}
  \item{rho}{Bandwidth parameter for the Gaussian kernel.}
  \item{degree}{Degree of the polynomial kernel.}
  \item{sigma2}{Variance parameter for the radial basis kernel.}
  \item{lambda}{Regularization parameter used during KQR fitting.}
}
\description{
Generates predictions for new data using a fitted Kernel Quantile Regression model. Computes the kernel between new inputs and training data, and applies the learned dual coefficients.
}
\details{
This function computes the prediction:
\deqn{
\hat{f}(x) = \frac{1}{\lambda} K(x, x_{\text{train}}) \theta,
}
where \eqn{K(x, x_{\text{train}})} is the kernel matrix between the new inputs and training data, and \eqn{\theta} is the vector of dual coefficients.
}
\value{
Returns a numeric vector of predicted values for each row in \code{x_new}.
}
\examples{
# Simulated training data
set.seed(123)
x_train <- matrix(sort(runif(20)), ncol = 1)
y <- sin(2 * pi * x_train) + rnorm(20, sd = 0.1)

# Define a simple Gaussian kernel function
gaussian_kernel <- function(x, y, rho = 10) {
  exp(-rho * sum((x - y)^2))
}

# Construct kernel matrix K (symmetric)
K <- matrix(0, nrow = 20, ncol = 20)
for (i in 1:20) {
  for (j in 1:20) {
    K[i, j] <- gaussian_kernel(x_train[i, ], x_train[j, ])
  }
}

# Fit theta using quantile LP (using custom fit function)
fit <- fit_kqr_lp(K, y, tau = 0.5, lambda = 1)

# Predict on a grid
x_grid <- matrix(seq(0, 1, length.out = 100), ncol = 1)
y_pred <- predict_kqr(x_new = x_grid, x_train = x_train, theta = fit$theta,
                      kernel = "gaussian", rho = 10, lambda = 1)

# Plot
plot(x_train, y)
lines(x_grid, y_pred, col = "blue", lwd = 2)
}
