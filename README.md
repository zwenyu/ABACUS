ABACUS
--------------------------------------------------

R package implementing change detection procedure in paper ABACUS: Unsupervised Multivariate Change Detection via Bayesian Source Separation (arXiv:1810.06167). ABACUS is a Bayesian source separation technique to recover latent signals while locating additive outliers and level shifts in the unsupervised setting.


Example
--------------------------------------------------

```r
set.seed(6)
dat = simData(P = 10, N = 100, s_s = 4, t0_s = 10, t1_s = 10, sps = TRUE)
# true additive outliers and level shifts
dat$times0
dat$times1
# true model components
dat$Y # observations matrix
dat$M # mixing matrix
t(dat$S) # latent matrix
dat$psi # noise variance

# run proposed method at K = 5
res = bsscpt(dat$Y, K = 5, burn_in = 500, numiter = 3000)
# estimated additive outliers and level shifts
res$cpt0
res$cpt1
# estimated model components
res$Mest
t(res$Sest)
res$psiest
```

Installation
--------------------------------------------------

``` r
# install.packages("devtools")
devtools::install_github("zwenyu/ABACUS")
```
