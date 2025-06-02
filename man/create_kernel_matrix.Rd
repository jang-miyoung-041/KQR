\name{create_kernel_mat}
\alias{create_kernel_mat}
\title{Create a Kernel Matrix}
\usage{
create_kernel_mat(x, kernel = "gaussian", rho = 1, degree = 3, sigma2 = 1)
}
\arguments{
  \item{x}{An \code{n x p} numeric matrix or a numeric vector of input data.}
  \item{kernel}{Type of kernel function to use. One of \code{"gaussian"}, \code{"poly"}, or \code{"rad"}.}
  \item{rho}{Bandwidth parameter for the Gaussian kernel.}
  \item{degree}{Degree of the polynomial kernel.}
  \item{sigma2}{Variance parameter for the radial basis kernel.}
}
\description{
Constructs an \code{n x n} kernel matrix based on the selected kernel function and input data.
}
\details{
This function computes the kernel matrix \eqn{K} where each element \eqn{K_{ij} = k(x_i, x_j)} for a given kernel function.
}
\value{
Returns an \code{n x n} symmetric kernel matrix.
}
\examples{
# Example with 5x2 matrix
set.seed(123)
x = matrix(rnorm(10), ncol = 2)

# Gaussian kernel
K1 = create_kernel_mat(x, kernel = "gaussian", rho = 0.5)

# Polynomial kernel
K2 = create_kernel_mat(x, kernel = "poly", degree = 2)

# Radial basis kernel
K3 = create_kernel_mat(x, kernel = "rad", sigma2 = 1)
}