\name{ABACUS-package}
\alias{ABACUS-package}
\alias{ABACUS}
\docType{package}
\title{\packageTitle{ABACUS}}
\description{\packageDescription{ABACUS}}
\author{
  \packageAuthor{ABACUS}
}
\references{
 	W. Zhang, D. Gilbert, D. Matteson. ABACUS: Unsupervised Multivariate Change Detection via Bayesian Source Separation. arXiv:1810.06167
}
\examples{
\dontrun{
set.seed(6)
dat = simData(P = 10, N = 100, s_s = 4, t0_s = 10, t1_s = 10, sps = TRUE)
# true additive outliers and level shifts
dat$times0
dat$times1

# run proposed method at K = 5
res = bsscpt(dat$Y, K = 5, burn_in = 500, numiter = 3000)
res$cpt0
res$cpt1
}
}