\name{fit_kqr_lp}
\alias{fit_kqr_lp}
\title{Fit Kernel Quantile Regression via Linear Programming (No Intercept)}
\usage{
fit_kqr_lp(K, y, tau = 0.5, lambda = 1)
}
\arguments{
  \item{K}{An \code{n x n} kernel matrix computed from training data.}
  \item{y}{A numeric response vector of length \code{n}.}
  \item{tau}{Quantile level between 0 and 1 (e.g., \code{0.5} for median regression).}
  \item{lambda}{Regularization parameter (positive scalar) controlling smoothness.}
}
\description{
Fits a Kernel Quantile Regression (KQR) model using linear programming to estimate the dual coefficients \eqn{\theta}, without including an intercept term.
}
\details{
This function solves the optimization problem:
\deqn{
\min_{\theta, \xi^+, \xi^-} \sum_{i=1}^n \tau \xi^+_i + (1 - \tau) \xi^-_i,
}
subject to:
\deqn{
y_i = (K \theta / \lambda)_i + \xi^+_i - \xi^-_i,
}
\deqn{
\theta_i \in [-(1 - \tau), \tau], \quad \xi^+_i, \xi^-_i \ge 0.
}
The kernel matrix \code{K} can be computed using a Gaussian, polynomial, or radial basis kernel.
}
\value{
A list with the following components:
\item{\code{theta}}{Estimated dual coefficients of length \code{n}.}
\item{\code{objval}}{Optimal value of the objective function.}
\item{\code{status}}{Solver status code. A value of \code{0} indicates successful convergence.}
}
\examples{
# Example usage with synthetic data
set.seed(1)
x <- matrix(runif(30), ncol = 1)
y <- sin(2 * pi * x) + rnorm(30, sd = 0.1)

# Kernel matrix (Gaussian)
gaussian_kernel <- function(x1, x2, rho = 10) {
  exp(-rho * sum((x1 - x2)^2))
}

K <- matrix(0, nrow = 30, ncol = 30)
for (i in 1:30) {
  for (j in 1:30) {
    K[i, j] <- gaussian_kernel(x[i, ], x[j, ])
  }
}

fit <- fit_kqr_lp(K, y, tau = 0.5, lambda = 1)
plot(x, y)
lines(x[order(x)], (K %*% fit$theta)[order(x)], col = "blue")
}
