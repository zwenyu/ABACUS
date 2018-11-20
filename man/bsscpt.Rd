\name{bsscpt}
\alias{bsscpt}
\title{
  bsscpt
}
\usage{
	bsscpt(Y, K = 5, burn_in = 500, numiter = 3000, a.psi = 1, b.psi = 1)
}
\arguments{
	\item{Y}{P-by-N observations matrix.}

	\item{K}{Upper bound of latent dimensionality.}

	\item{burn_in}{Burn-in for Gibbs sampling.}

	\item{numiter}{Number of iterations for Gibbs sampling.}

	\item{a.psi, b.psi}{Parameters for priors of noise variance.}
}
\description{
  \code{bsscpt} estimates components of the latent model and change point locations.
}