% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runBMCMC.power.full.R
\name{runBMCMC.power.full}
\alias{runBMCMC.power.full}
\title{Gaussian process with a full power-exponential (Full-Hybrid) correlation structure.}
\usage{
runBMCMC.power.full(nmcmc = 10000, burn = 4000, thin = 10, x, y, xtest,
  lambda.ini, lambda.w.ini, gamma.ini, gamma.w.ini)
}
\arguments{
\item{nmcmc}{Length of the MCMC chain. Default is 10,000.}

\item{burn}{Burn-in period. Default is 4,000.}

\item{thin}{Thinning parameter. Default is 10.}

\item{x}{Training set. Either a matrix or a data frame.}

\item{y}{True outputs of x.}

\item{xtest}{Hold-out set, of which the ouputs will be predicted. Must have the same dimension as x.}

\item{lambda.ini}{Initial values of the sensitivity parameters, at \eqn{\lambda} scale.}

\item{lambda.w.ini}{Initial values of the lengths (variances) for adaptive MCMC in estimating the sensitivity parameters, at \eqn{\lambda} scale.}

\item{gamma.ini}{Initial values of the smoothness parameters, at \eqn{\gamma} scale.}

\item{gamma.w.ini}{Initial values of the lengths (variances) for adaptive MCMC in estimating the smoothness parameters, at \eqn{\gamma} scale.}
}
\value{
pred.y return the predicted mean.

pred.var return the predicted variance.

accept.rate.lambda return a vector (1*d) of the acceptance rate of lambda.

accept.rate.gamma return a vector (1*d) of the acceptance rate of gamma.

mcmc.matrix.lambda return a matrix (nmcmc * d) of the sampled lambda.

mcmc.matrix.gamma return a matrix (nmcmc * d) of the sampled gamma.
}
\description{
The power-exponential structure is
\deqn{ R(x_i, x_j)= \exp(\sum_{k=1}^{d} \theta_k (x_{i, k} - x_{j,k})^{\alpha_k}),}
where \eqn{ 1 \le  i, j \le n}, \eqn{\theta_k > 0 } is the sensitivity parameter and
\eqn{ 1 \le \alpha_k \le 2 } is the smoothness parameter.
\deqn{ \lambda_k = \log_{e}(\rho_k/(1-\rho_k)),}
where \eqn{\rho_k = \exp(-\theta_k/4)}. After the above \eqn{\lambda} transformation, the sensitivity parameter
is unconstrained.
\deqn{ \gamma_k = -\exp(\frac{2-\alpha_k}{\alpha_k-1}).}
The smoothness parameter is also unconstrained at \eqn{\gamma} scale.
}
