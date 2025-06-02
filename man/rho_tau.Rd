\name{rho_tau}
\alias{rho_tau}
\title{Check Loss Function (Pinball Loss)}
\usage{
rho_tau(y, y_hat, tau = 0.6)
}
\arguments{
  \item{y}{A numeric vector of true values.}
  \item{y_hat}{A numeric vector of predicted values.}
  \item{tau}{A numeric value between 0 and 1 indicating the quantile level. Default is 0.6.}
}
\description{
Computes the check loss function, also known as the pinball loss, for quantile regression.
}
\details{
The check loss is defined as:
\deqn{
\rho_\tau(u) = \begin{cases}
\tau \cdot u & \text{if } u \ge 0 \\
(\tau - 1) \cdot u & \text{if } u < 0
\end{cases}
}
This function applies the formula element-wise to \code{y - y\_hat}.
}
\value{
A numeric vector of loss values, same length as \code{y}.
}
\examples{
# Example: compare predicted and true values at tau = 0.6
y = c(3.2, 4.5, 2.8)
y_hat = c(3.0, 4.8, 2.5)
rho_tau(y, y_hat, tau = 0.6)
}