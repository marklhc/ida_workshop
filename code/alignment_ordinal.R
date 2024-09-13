#m_items <- paste0("i", 1:5)
#dat <- readRDS("rds/dat.rds")

# fit a graded response model (equivalent to an ordinal cfa).
config_ord <- multipleGroup(data = dat[, m_items], group = dat$sample,
                            itemtype = "graded", invariance = "i3")
# extract coefficients
est_coef <- mirt.wrapper.coef(config_ord)$coef
# Remove coefficients for i3 for HSLS
est_coef[6, c("a1", "d1", "d2", "d3")] <- NA
# est_coef has columns group, item, a1, d1, d2, and d3. 
# a1 contains the loadings, d1-d3 are the thresholds. 
(els_coef <- est_coef[est_coef$group == "ELS",])
(hsls_coef <- est_coef[est_coef$group == "HSLS",])

# To compute threshold estimates for polytomous items, we’ll treat each
# threshold as an item in sirt::invariance.alignment so we modify the lambda and 
# weight matrices to be in compatible dimensions
 
# prepare lambda and nu matrices with estimates from the configural model
(lambda1 <- rbind(rep(els_coef$a1, 3), rep(hsls_coef$a1, 3)))
(nu1 <- rbind(unlist(els_coef[, c("d1", "d2", "d3")]), 
              unlist(hsls_coef[, c("d1", "d2", "d3")])))
# weight matrix
(wgt_mat <- as.matrix(rbind(
  colSums(!is.na(dat[dat$sample=="ELS", m_items])),
  colSums(!is.na(dat[dat$sample=="HSLS", m_items]))
)))

# perform alignment and obtain aligned parameters of the latent mean and latent 
# variance
ord_align <- invariance.alignment(lambda1, nu1, 
                                  wgt = sqrt(cbind(wgt_mat, wgt_mat, wgt_mat)))
ord_align$pars
ord_align$pars[2, 2]^2 # square the SD to get the variance, 0.9700555

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

# Make sure coefficients are comparable to aligned parameters
ord_align$lambda.aligned[, 1:5] - mirt.wrapper.coef(mod_aligned)$coef$a1
ord_align$nu.aligned - mirt.wrapper.coef(mod_aligned)$coef[c("d1", "d2", "d3")]

# Compute factor scores
fs_align_ord <- fscores(mod_aligned,
                        mean = c(0, 0),
                        cov = c(1, 1),
                        full.scores.SE = TRUE)

fs_align_ord <- cbind(fs_align_ord, dat[!all_na, "sample"])

saveRDS(fs_align_ord, 'rds/fs_align_ord.rds')


# As we observed covariances in the continuous case, we might want to explore if
# covariances may improve model fit.

# Estimate covariances between i1-i2, i2-i3, and i2-i4/capture the shared 
# variance between the items by estimating each as a factor (F2, F3, F4)
mirt_syn <- "
    F1 = i1, i2, i3, i4, i5
    F2 = i1, i2
    F3 = i2, i3
    F4 = i2, i4
    CONSTRAIN = (i1-i2, a2), (i2-i3, a3), (i2, i4, a4)
  "
mirtmodel_cov <- mirt.model(mirt_syn, itemnames = m_items)
config_ord2 <- multipleGroup(data = dat[, m_items], model = mirtmodel_cov, 
                             group = dat$sample, itemtype = "graded", 
                             invariance = "i3"
)

residuals(config_ord2, type = "Q3")  # identify strong unique covariances

# As the residual correlations are mostly negative, adding dimensions would 
# not help. We retain the analyses on the model without covariances. 
