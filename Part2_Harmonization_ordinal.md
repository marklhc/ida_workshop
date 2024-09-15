SREE Workshop Data Harmonization Illustrative Example: Ordinal Case
================
2024-09-18

# Install and load packages, prepare data.

``` r
# install and load frequently used packages
library(dplyr)
library(lavaan)
library(sirt)
library(mirt)
library(here)
# also install packages: haven, numDeriv
```

``` r
dir.create("rds", showWarnings = FALSE)
if (!file.exists("rds/dat.rds")) {
  source("code/download_data.R")
  source("code/prepare_data.R")
} else {
  dat <- readRDS("rds/dat.rds")
}
```

``` r
m_items <- paste0("i", 1:5)
# get subset of relevant variables
dat <- dat[, c("stu_id", "sample", "sex", "dropout", m_items)]
```

As a baseline for comparison against factor scores, we compute mean
scores. We opt for mean scores instead of sum scores due to the missing
item in HSLS.

``` r
dat$mean_score <- c(rowMeans(dat[dat$sample == "ELS", m_items], na.rm = TRUE),
                    rowMeans(dat[dat$sample == "HSLS", m_items[-3]], na.rm = TRUE))
```

# Obtain factor scores assuming ordinal data

## Step 1: Determine a partial invariance model through traditional MI modeling (ordinal case).

`traditional_mi_testing_ordinal.R` contains the specific steps we
followed to determine the final partial invariance model assuming
ordered categorical data. Five models were fit and compared to determine
a configural model followed by 24 additional models to determine the
final partial invariance model. For additional details, see Tse, W. W.
Y., Lai, M. H., & Zhang, Y. (2024). Does strict invariance matter? Valid
group mean comparisons with ordered-categorical items. Behavior Research
Methods, 56(4), 3117-3139.

Below is the final partial invariance model determined through this
procedure.

``` r
cfa_partial_ord <- 'group: ELS
                    math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5

                    i1 ~~ 1 * i1
                    i2 ~~ 1 * i2
                    i3 ~~ 1 * i3
                    i4 ~~ 1 * i4
                    i5 ~~ 1 * i5
                    i1 ~~ i2
                    i2 ~~ i3
                    i2 ~~ i4

                    math ~~ 1 * math
                    math ~ 0 * 1

                    i1 | th1 * t1
                    i1 | th2 * t2
                    i1 | th3 * t3

                    i2 | th4 * t1
                    i2 | th5 * t2
                    i2 | th6 * t3

                    i3 | th7 * t1
                    i3 | th8 * t2
                    i3 | th9 * t3

                    i4 | th10 * t1
                    i4 | th11 * t2
                    i4 | th12 * t3

                    i5 | th13 * t1
                    i5 | th14 * t2
                    i5 | th15 * t3

                    group: HSLS
                    math =~ NA * i1 + i2 + l4 * i4 + l5 * i5

                    i1 ~~ NA * i1
                    i2 ~~ NA * i2
                    i4 ~~ 1 * i4
                    i5 ~~ NA * i5
                    i1 ~~ i2
                    i2 ~~ i4

                    math ~~ NA * math
                    math ~ NA * 1

                    i1 | th1 * t1
                    i1 | th2a * t2
                    i1 | th3 * t3

                    i2 | th4 * t1
                    i2 | th5a * t2
                    i2 | th6 * t3

                    i4 | th10 * t1
                    i4 | th11a * t2
                    i4 | th12 * t3

                    i5 | th13 * t1
                    i5 | th14a * t2
                    i5 | th15a * t3
                    '

fit_partial_ord  <- cfa(cfa_partial_ord, data = dat, group = "sample",
                        estimator = "WLSMV", ordered = TRUE,
                        parameterization = "theta")
```

``` r
# extract cfa parameter estimates
est_partial_ord <- lavInspect(fit_partial_ord, what = "est")
```

## Step 2. Determine an approximate invariance model through alignment optimization (ordinal case).

``` r
# fit a graded response model (equivalent to an ordinal cfa).
config_ord <- multipleGroup(data = dat[, m_items], group = dat$sample,
                            itemtype = "graded", invariance = "i3")
```

    ## Warning: data contains response patterns with only NAs

    ## Iteration: 1, Log-Lik: -148893.119, Max-Change: 2.94493Iteration: 2, Log-Lik: -121681.876, Max-Change: 1.02598Iteration: 3, Log-Lik: -116632.176, Max-Change: 2.33449Iteration: 4, Log-Lik: -114397.320, Max-Change: 0.37731Iteration: 5, Log-Lik: -113564.107, Max-Change: 0.34777Iteration: 6, Log-Lik: -113213.615, Max-Change: 0.11836Iteration: 7, Log-Lik: -113019.312, Max-Change: 0.09583Iteration: 8, Log-Lik: -112890.480, Max-Change: 0.06939Iteration: 9, Log-Lik: -112798.159, Max-Change: 0.06192Iteration: 10, Log-Lik: -112727.795, Max-Change: 0.06161Iteration: 11, Log-Lik: -112669.044, Max-Change: 0.05656Iteration: 12, Log-Lik: -112628.204, Max-Change: 0.07860Iteration: 13, Log-Lik: -112592.222, Max-Change: 0.05452Iteration: 14, Log-Lik: -112563.804, Max-Change: 0.06469Iteration: 15, Log-Lik: -112540.323, Max-Change: 0.07189Iteration: 16, Log-Lik: -112520.386, Max-Change: 0.08542Iteration: 17, Log-Lik: -112503.335, Max-Change: 0.07462Iteration: 18, Log-Lik: -112488.102, Max-Change: 0.08043Iteration: 19, Log-Lik: -112475.719, Max-Change: 0.07494Iteration: 20, Log-Lik: -112464.134, Max-Change: 0.08256Iteration: 21, Log-Lik: -112454.230, Max-Change: 0.07295Iteration: 22, Log-Lik: -112444.860, Max-Change: 0.07793Iteration: 23, Log-Lik: -112436.878, Max-Change: 0.07163Iteration: 24, Log-Lik: -112429.242, Max-Change: 0.07180Iteration: 25, Log-Lik: -112422.770, Max-Change: 0.06431Iteration: 26, Log-Lik: -112416.313, Max-Change: 0.06590Iteration: 27, Log-Lik: -112411.012, Max-Change: 0.06551Iteration: 28, Log-Lik: -112405.418, Max-Change: 0.03419Iteration: 29, Log-Lik: -112398.136, Max-Change: 0.06105Iteration: 30, Log-Lik: -112390.050, Max-Change: 0.05553Iteration: 31, Log-Lik: -112383.794, Max-Change: 0.04784Iteration: 32, Log-Lik: -112381.207, Max-Change: 0.07588Iteration: 33, Log-Lik: -112376.713, Max-Change: 0.05402Iteration: 34, Log-Lik: -112375.687, Max-Change: 0.03883Iteration: 35, Log-Lik: -112373.311, Max-Change: 0.03430Iteration: 36, Log-Lik: -112370.115, Max-Change: 0.04521Iteration: 37, Log-Lik: -112367.277, Max-Change: 0.02007Iteration: 38, Log-Lik: -112365.142, Max-Change: 0.01484Iteration: 39, Log-Lik: -112363.503, Max-Change: 0.01506Iteration: 40, Log-Lik: -112356.328, Max-Change: 0.00971Iteration: 41, Log-Lik: -112355.728, Max-Change: 0.00845Iteration: 42, Log-Lik: -112355.268, Max-Change: 0.00733Iteration: 43, Log-Lik: -112353.351, Max-Change: 0.00856Iteration: 44, Log-Lik: -112353.198, Max-Change: 0.00411Iteration: 45, Log-Lik: -112353.070, Max-Change: 0.00383Iteration: 46, Log-Lik: -112352.533, Max-Change: 0.00276Iteration: 47, Log-Lik: -112352.488, Max-Change: 0.00219Iteration: 48, Log-Lik: -112352.454, Max-Change: 0.00208Iteration: 49, Log-Lik: -112352.317, Max-Change: 0.00138Iteration: 50, Log-Lik: -112352.300, Max-Change: 0.00068Iteration: 51, Log-Lik: -112352.295, Max-Change: 0.00115Iteration: 52, Log-Lik: -112352.278, Max-Change: 0.00116Iteration: 53, Log-Lik: -112352.271, Max-Change: 0.00102Iteration: 54, Log-Lik: -112352.265, Max-Change: 0.00030Iteration: 55, Log-Lik: -112352.265, Max-Change: 0.00024Iteration: 56, Log-Lik: -112352.264, Max-Change: 0.00087Iteration: 57, Log-Lik: -112352.261, Max-Change: 0.00067Iteration: 58, Log-Lik: -112352.255, Max-Change: 0.00029Iteration: 59, Log-Lik: -112352.254, Max-Change: 0.00058Iteration: 60, Log-Lik: -112352.252, Max-Change: 0.00014Iteration: 61, Log-Lik: -112352.252, Max-Change: 0.00068Iteration: 62, Log-Lik: -112352.251, Max-Change: 0.00034Iteration: 63, Log-Lik: -112352.250, Max-Change: 0.00048Iteration: 64, Log-Lik: -112352.248, Max-Change: 0.00019Iteration: 65, Log-Lik: -112352.248, Max-Change: 0.00049Iteration: 66, Log-Lik: -112352.247, Max-Change: 0.00012Iteration: 67, Log-Lik: -112352.247, Max-Change: 0.00050Iteration: 68, Log-Lik: -112352.246, Max-Change: 0.00022Iteration: 69, Log-Lik: -112352.245, Max-Change: 0.00043Iteration: 70, Log-Lik: -112352.245, Max-Change: 0.00024Iteration: 71, Log-Lik: -112352.244, Max-Change: 0.00044Iteration: 72, Log-Lik: -112352.244, Max-Change: 0.00019Iteration: 73, Log-Lik: -112352.243, Max-Change: 0.00016Iteration: 74, Log-Lik: -112352.243, Max-Change: 0.00036Iteration: 75, Log-Lik: -112352.243, Max-Change: 0.00061Iteration: 76, Log-Lik: -112352.242, Max-Change: 0.00030Iteration: 77, Log-Lik: -112352.242, Max-Change: 0.00056Iteration: 78, Log-Lik: -112352.242, Max-Change: 0.00024Iteration: 79, Log-Lik: -112352.242, Max-Change: 0.00020Iteration: 80, Log-Lik: -112352.241, Max-Change: 0.00037Iteration: 81, Log-Lik: -112352.241, Max-Change: 0.00016Iteration: 82, Log-Lik: -112352.241, Max-Change: 0.00013Iteration: 83, Log-Lik: -112352.241, Max-Change: 0.00026Iteration: 84, Log-Lik: -112352.241, Max-Change: 0.00010Iteration: 85, Log-Lik: -112352.241, Max-Change: 0.00044Iteration: 86, Log-Lik: -112352.241, Max-Change: 0.00018Iteration: 87, Log-Lik: -112352.240, Max-Change: 0.00035Iteration: 88, Log-Lik: -112352.240, Max-Change: 0.00023Iteration: 89, Log-Lik: -112352.240, Max-Change: 0.00009

``` r
# extract coefficients
est_coef <- mirt.wrapper.coef(config_ord)$coef
# Remove coefficients for i3 for HSLS
est_coef[6, c("a1", "d1", "d2", "d3")] <- NA
# est_coef has columns group, item, a1, d1, d2, and d3.
# a1 contains the loadings, d1-d3 are the thresholds.
(els_coef <- est_coef[est_coef$group == "ELS",])
```

    ##   group item       a1       d1         d2        d3
    ## 1   ELS   i1 3.367244 4.716831 -0.4397776 -3.074938
    ## 3   ELS   i2 3.483662 3.490577 -0.8780415 -4.181195
    ## 5   ELS   i3 4.365504 4.500842 -0.4857222 -4.260119
    ## 7   ELS   i4 4.372924 5.632033  0.2849476 -3.562358
    ## 9   ELS   i5 4.337701 5.732986  0.4396606 -3.412799

``` r
(hsls_coef <- est_coef[est_coef$group == "HSLS",])
```

    ##    group item       a1       d1       d2        d3
    ## 2   HSLS   i1 3.863295 7.498264 3.094064 -2.864859
    ## 4   HSLS   i2 2.996175 4.982647 1.199976 -3.369710
    ## 6   HSLS   i3       NA       NA       NA        NA
    ## 8   HSLS   i4 4.313544 8.691290 4.301483 -2.717268
    ## 10  HSLS   i5 3.824915 7.657129 3.279324 -3.066883

``` r
# To compute threshold estimates for polytomous items, we’ll treat each
# threshold as an item in sirt::invariance.alignment so we modify the lambda and
# weight matrices to be in compatible dimensions

# prepare lambda and nu matrices with estimates from the configural model
(lambda1 <- rbind(rep(els_coef$a1, 3), rep(hsls_coef$a1, 3)))
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]     [,9]
    ## [1,] 3.367244 3.483662 4.365504 4.372924 4.337701 3.367244 3.483662 4.365504 4.372924
    ## [2,] 3.863295 2.996175       NA 4.313544 3.824915 3.863295 2.996175       NA 4.313544
    ##         [,10]    [,11]    [,12]    [,13]    [,14]    [,15]
    ## [1,] 4.337701 3.367244 3.483662 4.365504 4.372924 4.337701
    ## [2,] 3.824915 3.863295 2.996175       NA 4.313544 3.824915

``` r
(nu1 <- rbind(unlist(els_coef[, c("d1", "d2", "d3")]),
              unlist(hsls_coef[, c("d1", "d2", "d3")])))
```

    ##           d11      d12      d13      d14      d15        d21        d22        d23
    ## [1,] 4.716831 3.490577 4.500842 5.632033 5.732986 -0.4397776 -0.8780415 -0.4857222
    ## [2,] 7.498264 4.982647       NA 8.691290 7.657129  3.0940640  1.1999761         NA
    ##            d24       d25       d31       d32       d33       d34       d35
    ## [1,] 0.2849476 0.4396606 -3.074938 -4.181195 -4.260119 -3.562358 -3.412799
    ## [2,] 4.3014829 3.2793241 -2.864859 -3.369710        NA -2.717268 -3.066883

``` r
# weight matrix
(wgt_mat <- as.matrix(rbind(
  colSums(!is.na(dat[dat$sample=="ELS", m_items])),
  colSums(!is.na(dat[dat$sample=="HSLS", m_items]))
)))
```

    ##         i1    i2    i3    i4    i5
    ## [1,] 11391 11432 11047 10823 10661
    ## [2,] 19049 19006     0 18926 18976

``` r
# perform alignment and obtain aligned parameters of the latent mean and latent
# variance
ord_align <- invariance.alignment(lambda1, nu1,
                                  wgt = sqrt(cbind(wgt_mat, wgt_mat, wgt_mat)))
ord_align$pars
```

    ##      alpha0     psi0
    ## G1 0.000000 1.000000
    ## G2 0.699614 0.984914

``` r
ord_align$pars[2, 2]^2 # square the SD to get the variance, 0.9700555
```

    ## [1] 0.9700555

``` r
# Constrain latent means and variances to the aligned levels
mirtmodel_al <- "
    F1 = i1, i2, i3, i4, i5
    START [HSLS] = (GROUP, MEAN_1, 0.699614), (GROUP, COV_11, 0.9700555)
  "
# Specify the mirt model explicitly providing the item names
mirtmodel_al <- mirt.model(mirtmodel_al, itemnames = m_items)
# Drop observations where all math items are NAs to avoid issues with `fscores()`
all_na <- rowSums(!is.na(dat[, m_items])) == 0

# Fit the aligned model
mod_aligned <- multipleGroup(dat[!all_na, m_items],
                             model = mirtmodel_al,
                             group = dat$sample[!all_na],
                             itemtype = "graded",
                             invariance = "i3")
```

    ## Iteration: 1, Log-Lik: -149122.928, Max-Change: 2.99206Iteration: 2, Log-Lik: -122268.510, Max-Change: 1.00989Iteration: 3, Log-Lik: -116604.720, Max-Change: 0.88186Iteration: 4, Log-Lik: -114771.156, Max-Change: 0.56024Iteration: 5, Log-Lik: -113998.474, Max-Change: 0.34573Iteration: 6, Log-Lik: -113576.493, Max-Change: 0.14773Iteration: 7, Log-Lik: -113338.177, Max-Change: 0.11380Iteration: 8, Log-Lik: -113164.620, Max-Change: 0.11527Iteration: 9, Log-Lik: -113028.697, Max-Change: 0.09718Iteration: 10, Log-Lik: -112925.143, Max-Change: 0.10747Iteration: 11, Log-Lik: -112840.641, Max-Change: 0.13377Iteration: 12, Log-Lik: -112761.400, Max-Change: 0.09591Iteration: 13, Log-Lik: -112712.165, Max-Change: 0.08627Iteration: 14, Log-Lik: -112668.741, Max-Change: 0.06478Iteration: 15, Log-Lik: -112635.475, Max-Change: 0.08256Iteration: 16, Log-Lik: -112606.966, Max-Change: 0.06835Iteration: 17, Log-Lik: -112578.995, Max-Change: 0.10278Iteration: 18, Log-Lik: -112556.879, Max-Change: 0.06556Iteration: 19, Log-Lik: -112535.560, Max-Change: 0.06207Iteration: 20, Log-Lik: -112518.247, Max-Change: 0.08714Iteration: 21, Log-Lik: -112501.625, Max-Change: 0.07367Iteration: 22, Log-Lik: -112488.748, Max-Change: 0.09351Iteration: 23, Log-Lik: -112476.673, Max-Change: 0.04627Iteration: 24, Log-Lik: -112466.929, Max-Change: 0.07436Iteration: 25, Log-Lik: -112455.428, Max-Change: 0.04324Iteration: 26, Log-Lik: -112446.584, Max-Change: 0.05546Iteration: 27, Log-Lik: -112437.554, Max-Change: 0.03356Iteration: 28, Log-Lik: -112428.673, Max-Change: 0.04573Iteration: 29, Log-Lik: -112420.469, Max-Change: 0.03291Iteration: 30, Log-Lik: -112410.900, Max-Change: 0.04039Iteration: 31, Log-Lik: -112405.303, Max-Change: 0.02947Iteration: 32, Log-Lik: -112399.990, Max-Change: 0.03447Iteration: 33, Log-Lik: -112395.274, Max-Change: 0.02799Iteration: 34, Log-Lik: -112391.043, Max-Change: 0.03011Iteration: 35, Log-Lik: -112387.298, Max-Change: 0.02734Iteration: 36, Log-Lik: -112383.858, Max-Change: 0.02437Iteration: 37, Log-Lik: -112374.317, Max-Change: 0.02085Iteration: 38, Log-Lik: -112372.361, Max-Change: 0.01655Iteration: 39, Log-Lik: -112369.886, Max-Change: 0.01497Iteration: 40, Log-Lik: -112359.136, Max-Change: 0.01205Iteration: 41, Log-Lik: -112358.159, Max-Change: 0.00980Iteration: 42, Log-Lik: -112357.396, Max-Change: 0.00886Iteration: 43, Log-Lik: -112354.194, Max-Change: 0.00613Iteration: 44, Log-Lik: -112353.885, Max-Change: 0.00525Iteration: 45, Log-Lik: -112353.670, Max-Change: 0.00473Iteration: 46, Log-Lik: -112352.917, Max-Change: 0.00434Iteration: 47, Log-Lik: -112352.793, Max-Change: 0.00349Iteration: 48, Log-Lik: -112352.716, Max-Change: 0.00313Iteration: 49, Log-Lik: -112352.431, Max-Change: 0.00329Iteration: 50, Log-Lik: -112352.387, Max-Change: 0.00116Iteration: 51, Log-Lik: -112352.373, Max-Change: 0.00193Iteration: 52, Log-Lik: -112352.330, Max-Change: 0.00095Iteration: 53, Log-Lik: -112352.321, Max-Change: 0.00088Iteration: 54, Log-Lik: -112352.313, Max-Change: 0.00182Iteration: 55, Log-Lik: -112352.290, Max-Change: 0.00201Iteration: 56, Log-Lik: -112352.284, Max-Change: 0.00078Iteration: 57, Log-Lik: -112352.278, Max-Change: 0.00054Iteration: 58, Log-Lik: -112352.273, Max-Change: 0.00036Iteration: 59, Log-Lik: -112352.272, Max-Change: 0.00084Iteration: 60, Log-Lik: -112352.269, Max-Change: 0.00091Iteration: 61, Log-Lik: -112352.263, Max-Change: 0.00053Iteration: 62, Log-Lik: -112352.259, Max-Change: 0.00037Iteration: 63, Log-Lik: -112352.258, Max-Change: 0.00022Iteration: 64, Log-Lik: -112352.258, Max-Change: 0.00059Iteration: 65, Log-Lik: -112352.255, Max-Change: 0.00057Iteration: 66, Log-Lik: -112352.253, Max-Change: 0.00016Iteration: 67, Log-Lik: -112352.253, Max-Change: 0.00054Iteration: 68, Log-Lik: -112352.252, Max-Change: 0.00066Iteration: 69, Log-Lik: -112352.251, Max-Change: 0.00030Iteration: 70, Log-Lik: -112352.251, Max-Change: 0.00011Iteration: 71, Log-Lik: -112352.250, Max-Change: 0.00050Iteration: 72, Log-Lik: -112352.249, Max-Change: 0.00032Iteration: 73, Log-Lik: -112352.249, Max-Change: 0.00017Iteration: 74, Log-Lik: -112352.248, Max-Change: 0.00046Iteration: 75, Log-Lik: -112352.248, Max-Change: 0.00048Iteration: 76, Log-Lik: -112352.247, Max-Change: 0.00020Iteration: 77, Log-Lik: -112352.247, Max-Change: 0.00009

``` r
est_align_ord <- mirt.wrapper.coef(mod_aligned)$coef
est_align_ord[6, c("a1", "d1", "d2", "d3")] <- NA
```

``` r
# Make sure coefficients are comparable to aligned parameters
ord_align$lambda.aligned[, 1:5] - est_align_ord$a1
```

    ##               I1           I2           I3            I4            I5
    ## G1 -6.425054e-05 2.104800e-04 9.607714e-05 -0.0006018592 -6.585290e-04
    ## G2  8.397124e-06 5.518844e-05           NA  0.0000146369  2.559408e-05

``` r
est_align_ord_th <- matrix(unlist(est_align_ord[, c("d1", "d2", "d3")]),
                           nrow = 2)
est_align_ord_lam <- matrix(est_align_ord[, c("a1")], nrow = 2)
```

## Step 3: Compute factor scores and obtain standard errors

In the ordinal case, we compute scores from the partial invariance model
using the EBM (Empirical Bayes Modal) approach, which is one of the two
options available in `lavaan::lavPredict()` for categorical data. We
compute EAP (expected a-posteriori) factor scores from the approximate
invariance model, which is the default option for ’mirt::fscores()\`

``` r
fs_partial_ord <- lavPredict(fit_partial_ord, method = "EBM", se = TRUE)
```

    ## Warning in lav_predict_internal(lavmodel = lavmodel, lavdata = lavdata, : lavaan
    ## WARNING: standard errors not available (yet) for non-normal data

``` r
fs_align_ord <- as.data.frame(fscores(mod_aligned, mean = c(0, 0),
                                      cov = c(1, 1), full.scores.SE = TRUE,
                                      method = "EAP"))
# store the factor scores
score_df_ord <- cbind(dat, approx_ord = NA, approx_ord_SE = NA)
score_df_ord$approx_ord[!all_na] <- fs_align_ord[, "F1"]
score_df_ord$approx_ord_SE[!all_na] <- fs_align_ord[, "SE_F1"]
```

Note that SE are not available for non-normal data in lavaan so we
cannot obtain SE for the factor scores computed using the partial
invariance model assuming ordinal data.

## Step 4: Compute reliability.

EAP score reliability
$$\text{EAP score reliability}=1- \frac{\text{SE}^2}{\psi}$$

``` r
# obtain latent variances for the two groups
psi_align_ord_ELS <- mirt.wrapper.coef(mod_aligned)$GroupPars$ELS[2]
psi_align_ord_HSLS <- mirt.wrapper.coef(mod_aligned)$GroupPars$HSLS[2]

rel_approx_ord_ELS <- 1 - score_df_ord[score_df_ord$sample == "ELS",
                                       "approx_ord_SE"]^2 /  psi_align_ord_ELS

rel_approx_ord_HSLS <- 1 - score_df_ord[score_df_ord$sample == "HSLS",
                                        "approx_ord_SE"]^2 / psi_align_ord_HSLS
score_df_ord$approx_ord_rel <- c(rel_approx_ord_ELS, rel_approx_ord_HSLS)
score_df_ord[,3:6] <- apply(score_df_ord[,3:6], FUN = as.numeric, MARGIN = 2)
```

``` r
# store error variances in a data frame
score_df_ord$approx_ord_ev <- score_df_ord$approx_ord_SE^2 * score_df_ord$approx_ord_rel
head(score_df_ord, 2)
```

    ##   stu_id sample sex dropout i1 i2 i3 i4 i5 mean_score approx_ord approx_ord_SE
    ## 1 101101    ELS   0       0  2  1  2  2  1        1.6  -1.018915     0.2476258
    ## 2 101102    ELS   0       0  4  3  4  4  4        3.8   1.312507     0.2987173
    ##   approx_ord_rel approx_ord_ev
    ## 1      0.9386815    0.05755857
    ## 2      0.9107680    0.08126969

``` r
saveRDS(score_df_ord, "rds/score_df_ord.rds")
saveRDS(est_partial_ord, "rds/est_partial_ord.rds")
saveRDS(est_align_ord, "rds/est_align_ord.rds")
saveRDS(est_align_ord_th, "rds/est_align_ord_th.rds")
saveRDS(est_align_ord_lam, "rds/est_align_ord_lam.rds")
saveRDS(mod_aligned, "rds/mod_aligned.rds")
saveRDS(fs_partial_ord, "rds/fs_partial_ord.rds")
```
