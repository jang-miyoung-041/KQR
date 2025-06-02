\name{predict_kqr}
\alias{predict_kqr}
\title{Prediction for Kernel Quantile Regression}
\usage{
predict_kqr(alpha, x_train, x_new, kernel = "gaussian", param = list(rho = 1))
}
\arguments{
  \item{alpha}{Numeric vector of fitted dual coefficients (usually obtained from kernel quantile regression).}
  \item{x_train}{Matrix of training data used to estimate \code{alpha}, with observations in rows.}
  \item{x_new}{Matrix of new input data points to predict, with observations in rows.}
  \item{kernel}{Kernel type to use: \code{"gaussian"}, \code{"poly"}, or \code{"rbf"}.}
  \item{param}{A list of kernel parameters. For example, \code{list(rho = 1)} for Gaussian kernel, \code{list(degree = 2)} for polynomial, or \code{list(sigma2 = 1)} for RBF kernel.}
}
\description{
Generates predictions for new input data using kernel quantile regression and pre-estimated dual coefficients.
}
\details{
The function computes the predicted value at a new input \code{x} using the formula:
\deqn{
f(x) = \sum_{i=1}^n \alpha_i K(x, x_i)
}
where \eqn{K} is the specified kernel function.
}
\value{
A numeric vector of predicted values for each row of \code{x_new}.
}
\examples{
set.seed(1)
x_train = matrix(rnorm(10), ncol = 1)
x_new   = matrix(seq(-2, 2, length.out = 5), ncol = 1)
alpha   = rnorm(nrow(x_train))

# Gaussian kernel prediction
predict_kqr(alpha, x_train, x_new, kernel = "gaussian", param = list(rho = 0.5))
}
