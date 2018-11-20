\name{prettycpt}
\alias{prettycpt}
\title{
  prettycpt
}
\usage{
	prettycpt(S, cpt0, cpt1)
}
\arguments{
	\item{S}{Estimated K-by-N latent matrix.}

	\item{cpt0}{Estimated locations of additive outliers.}

	\item{cpt1}{Estimated locations of level shifts.}
}
\description{
  \code{prettycpt} uses dynamic programming to prune the change points if they are overestimated due to data violation of model assumptions. Pruned sets of additive outliers and level shifts are returned.
}