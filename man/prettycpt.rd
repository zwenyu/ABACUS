\name{prettycpt}
\alias{prettycpt}
\title{
  prettycpt
}
\usage{
	prettycpt(S, cpt0, cpt1)
}
\arguments{
	\item{S}{Estimated latent matrix.}

	\item{cpt0}{Estimated locations of additive outliers.}

	\item{cpt1}{Estimated locations of level shifts.}
}
\description{
  \code{prettycpt} uses dynamic programming to prune the set of change points, if changes are overestimated due to data violation of model assumptions.
}