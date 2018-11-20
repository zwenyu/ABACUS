\name{simData}
\alias{simData}
\title{
  simData
}
\usage{
	simData(P, N, s_s, t0_s, t1_s, sps = FALSE, times = NULL)
}
\arguments{
	\item{P}{Dimensionality of observations.}

	\item{N}{Sample size.}

	\item{s_s}{Number of latent signals with changes.}

	\item{t0_s}{Number of additive outliers.}

	\item{t1_s}{Number of level shifts.}

	\item{sps}{Logical. If \code{TRUE}, not all latent signals change together.}
	
	\item{times}{Change locations. If \code{NULL}, locations are sampled uniformly at random.}
}
\description{
  \code{simData} simulates a P-by-N observations matrix with the specified properties. Model components and change locations are also returned.
}