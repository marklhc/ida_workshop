# repeat analyses for d_noNAs dataset

# assuming libraries have been loaded and the data files have been read in
# (note that the files were created by running `Part1_Harmonization_continuous.Rmd`
# and `Part1_Harmonization_ordinal.Rmd` so if these files have not been run please
# run them first.
dat <- readRDS("rds/dat.rds")
m_items <- paste0("i", 1:5)
complete_cases <- complete.cases(dat[, m_items]) & dat$sample == "ELS" |
  complete.cases(dat[, m_items[-3]]) & dat$sample == "HSLS"
d_noNAs <- dat[complete_cases, ]

###### Mean scores (M_B2.)
d_noNAs$mean_score <- 
  c(rowMeans(d_noNAs[d_noNAs$sample == "ELS", m_items], na.rm = TRUE),
    rowMeans(d_noNAs[d_noNAs$sample == "HSLS", m_items[-3]], na.rm = TRUE))

dim(d_noNAs)
score_df_noNAs <- as.data.frame(cbind(
  "sample" = d_noNAs[,"sample"], "dropout" = d_noNAs$dropout,
  "mean_score" = d_noNAs$mean_score[!is.na(d_noNAs$mean_score)]))
head(score_df_noNAs, 2)




###### Partial invariance (M4)

cfa_partial <-   'group: ELS
                  math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                         
                  i1 ~ NA * 1
                  i2 ~ NA * 1
                  i3 ~ nu3 * 1
                  i4 ~ nu4 * 1
                  i5 ~ nu5 * 1
                         
                  i1 ~~ theta1_1 * i1
                  i2 ~~ theta2_1 * i2
                  i3 ~~ theta3 * i3
                  i4 ~~ theta4_1 * i4
                  i5 ~~ theta5_1 * i5
                  i1 ~~ i2
                  i2 ~~ cov3 * i3
                  i2 ~~ i4
                         
                  math ~~ 1 * math
                  math ~ 0 * 1      
                     
                  group: HSLS
                  math =~ NA * i1 + i2 + l4 * i4 + l5 * i5
                         
                  i1 ~ NA * 1
                  i2 ~ NA * 1
                  i4 ~ nu4 * 1
                  i5 ~ nu5 * 1
                         
                  i1 ~~ theta1_2 * i1
                  i2 ~~ theta2_2 * i2
                  i4 ~~ theta4_2 * i4
                  i5 ~~ theta5_2 * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  '
if (!file.exists("rds/fit_partial_noNAs.rds")) {
  # Note that MI testing with dataset d_noNAs led to the same partial model and 
  # we do not provide the code for the procedure separately for brevity. 
  # If interested, please substitute 'd_noNAs' for 'dat' as the dataset and
  # rerun `traditional_mi_testing.R`.
  fit_partial_noNAs <- cfa(cfa_partial, data = d_noNAs, group = "sample", 
                           estimator = "MLR", missing = "FIML", se = "robust.mlr")
  est_partial_noNAs  <- lavInspect(fit_partial_noNAs , what = "est")
  saveRDS(fit_partial_noNAs, "rds/fit_partial_noNAs.rds")
  saveRDS(est_partial_noNAs, "rds/est_partial_noNAs.rds")
} else {
  fit_partial_noNAs <- readRDS("rds/fit_partial_noNAs.rds")
  est_partial_noNAs <- readRDS("rds/est_partial_noNAs.rds")
}

###### Alignment optimization (M5)
cfa_config <-  'group: ELS
                math =~ NA * i1 + l2_1 * i2 + l3 * i3 + l4_1 * i4 + l5_1 * i5
                       
                i1 ~ nu1_1 * 1
                i2 ~ nu2_1 * 1
                i3 ~ nu3 * 1
                i4 ~ nu4_1 * 1
                i5 ~ nu5_1 * 1
                       
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2_1 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4_1 * i4
                i5 ~~ theta5_1 * i5
                i1 ~~ i2
                i2 ~~ cov3 * i3
                i2 ~~ i4
                       
                math ~~ 1 * math
                math ~ 0 * 1      
                   
                group: HSLS
                math =~ NA * i1 + l2_2 * i2 + l4_2 * i4 + l5_2 * i5
                       
                i1 ~ nu1_2 * 1
                i2 ~ nu2_2 * 1
                i4 ~ nu4_2 * 1
                i5 ~ nu5_2 * 1
                       
                i1 ~~ theta1_2 * i1
                i2 ~~ theta2_2 * i2
                i4 ~~ theta4_2 * i4
                i5 ~~ theta5_2 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ 1 * math
                math ~ 0 * 1  
                '
fit_config_noNAs  <- cfa(cfa_config, data = d_noNAs, group = "sample", 
                         estimator = "MLR", missing = "FIML", se = "robust.mlr")
# extract the measurement intercept and loading estimates from the configural model
est_config_noNAs <- lavInspect(fit_config_noNAs, what = "est")
config_lambda_noNAs <- merge(est_config_noNAs$ELS$lambda, 
                             est_config_noNAs$HSLS$lambda,
                             by = "row.names", all = TRUE)
config_lambda_mat_noNAs <- t(config_lambda_noNAs[, -1])
config_nu_noNAs <- merge(est_config_noNAs$ELS$nu, est_config_noNAs$HSLS$nu,
                         by = "row.names", all = TRUE)
config_nu_mat_noNAs <- t(config_nu_noNAs[, -1])

# obtain a weight matrix of the sample sizes for each item by sample
(wgt_mat_noNAs <- as.matrix(rbind(
  colSums(!is.na(d_noNAs[d_noNAs$sample=="ELS", m_items])),
  colSums(!is.na(d_noNAs[d_noNAs$sample=="HSLS", m_items]))
)))

aligned_par_noNAs <- invariance.alignment(lambda = config_lambda_mat_noNAs,
                                          nu = config_nu_mat_noNAs,
                                          fixed = TRUE, # fix SD of first group to one
                                          wgt = sqrt(wgt_mat_noNAs))

# modify the configural model specification string to freely estimate the latent
# mean and variance in both groups, and substitute in the aligned estimates from
# the configural model.
cfa_align_noNAs <-  
               'group: ELS
                math =~ l1_1 * i1 + l2_1 * i2 + l3 * i3 + l4_1 * i4 + l5_1 * i5
                       
                i1 ~ nu1_1 * 1
                i2 ~ nu2_1 * 1
                i3 ~ nu3 * 1
                i4 ~ nu4_1 * 1
                i5 ~ nu5_1 * 1
                       
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2_1 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4_1 * i4
                i5 ~~ theta5_1 * i5
                i1 ~~ i2
                i2 ~~ cov3 * i3
                i2 ~~ i4

                math ~~ NA * math
                math ~ NA * 1
                
                group: HSLS
                math =~ l1_2 * i1 + l2_2 * i2 + l4_2 * i4 + l5_2 * i5
                       
                i1 ~ nu1_2 * 1
                i2 ~ nu2_2 * 1
                i4 ~ nu4_2 * 1
                i5 ~ nu5_2 * 1
                       
                i1 ~~ theta1_2 * i1
                i2 ~~ theta2_2 * i2
                i4 ~~ theta4_2 * i4
                i5 ~~ theta5_2 * i5
                i1 ~~ i2
                i2 ~~ i4
                
                math ~~ NA * math
                math ~ NA * 1
                '
cfa_align_noNAs <- sub("l1_1", paste0(aligned_par_noNAs$lambda.aligned[,1][1], 
                                      collapse = ","), cfa_align_noNAs)
cfa_align_noNAs <- sub("l1_2", paste0(aligned_par_noNAs$lambda.aligned[,1][2], 
                                      collapse = ","), cfa_align_noNAs)
cfa_align_noNAs <- gsub("nu1_1", paste0(aligned_par_noNAs$nu.aligned[,1][1],
                                        collapse = ","), cfa_align_noNAs)
cfa_align_noNAs <- gsub("nu1_2", paste0(aligned_par_noNAs$nu.aligned[,1][2],
                                        collapse = ","), cfa_align_noNAs)

if (!file.exists("rds/fit_align_noNAs.rds")| !file.exists("rds/est_align_noNAs.rds")) {
  fit_align_noNAs  <- cfa(model = cfa_align_noNAs,  
                          data = d_noNAs,  estimator = "MLR", group = "sample",
                          missing = "FIML", se = "robust.mlr")
  est_align_noNAs <- lavInspect(fit_align_noNAs, what = "est")
  saveRDS(fit_align_noNAs, "rds/fit_align_noNAs.rds")
  saveRDS(est_align_noNAs, "rds/est_align_noNAs.rds")
} else {
  fit_align_noNAs <- readRDS("rds/fit_align_noNAs.rds")
  est_align_noNAs <- readRDS("rds/est_align_noNAs.rds")
}

# compute factor scores for the continuous case
fs_partial_noNAs <- lavPredict(fit_partial_noNAs, method = "bartlett", se = TRUE)
fs_align_noNAs <- lavPredict(fit_align_noNAs, method = "bartlett", se = TRUE)

# obtain standard errors of factor scores for each group
fs_partial_noNAs_SE <- unlist(attributes(fs_partial_noNAs)$se)
fs_align_noNAs_SE <- unlist(attributes(fs_align_noNAs)$se)

# store FS and FS SE in the score_df_noNAs data frame
score_df_noNAs <- as.data.frame(cbind(score_df_noNAs,
  "partial_cont" = as.numeric(unlist(fs_partial_noNAs)[!is.na(unlist(fs_partial_noNAs))]),
  "partial_SE" = fs_partial_noNAs_SE[!is.na(fs_partial_noNAs_SE)],
  "approx_cont" = as.numeric(unlist(fs_align_noNAs)[!is.na(unlist(fs_align_noNAs))]),
  "approx_SE" = fs_align_noNAs_SE[!is.na(fs_align_noNAs_SE)]))                                  

# compute reliability
bartlett_rel <- function(psi, SE) { psi / (psi + SE^2) } 

# compute reliability for each individual
rel_partial_noNAs_ELS <- 
  bartlett_rel(c(est_partial_noNAs$ELS$psi), 
               score_df_noNAs[score_df_noNAs$sample == "ELS", "partial_SE"])
rel_partial_noNAs_HSLS <- 
  bartlett_rel(c(est_partial_noNAs$HSLS$psi), 
               score_df_noNAs[score_df_noNAs$sample == "HSLS", "partial_SE"])
score_df_noNAs$partial_rel <- c(rel_partial_noNAs_ELS, rel_partial_noNAs_HSLS)

rel_approx_noNAs_ELS <- 
  bartlett_rel(c(est_align_noNAs$ELS$psi), 
               score_df_noNAs[score_df_noNAs$sample == "ELS", "approx_SE"])
rel_approx_noNAs_HSLS <- 
  bartlett_rel(c(est_align_noNAs$HSLS$psi),
               score_df_noNAs[score_df_noNAs$sample == "HSLS", "approx_SE"])
score_df_noNAs$approx_rel <- c(rel_approx_noNAs_ELS, rel_approx_noNAs_HSLS)

# store error variances
score_df_noNAs <- as.data.frame(cbind(score_df_noNAs,
  "partial_ev" = score_df_noNAs$partial_SE^2,
  "approx_ev" = score_df_noNAs$approx_SE))  

# extract latent variances
# psi_partial_ELS_noNAs <- c(est_partial_noNAs$ELS$psi)
# psi_partial_HSLS_noNAs <- c(est_partial_noNAs$HSLS$psi)
# psi_align_ELS_noNAs <- c(est_align_noNAs$ELS$psi)
# psi_align_HSLS_noNAs <- c(est_align_noNAs$HSLS$psi)

####### MI testing for the ordinal case (M6) 
# By default the rows with missing data are dropped.
# See `traditional_mi_testing_ordinal.R` for details

# read in the factor scores computed in `Part2_Harmonization_ordinal.Rmd`
fs_partial_ord <- readRDS("rds/fs_partial_ord.rds")

####### Alignment optimization (ordinal case; M7)
# fit a graded response model (equivalent to an ordinal cfa).
config_ord_noNAs <- multipleGroup(data = d_noNAs[, m_items], 
                                  group = d_noNAs$sample,
                                  itemtype = "graded", invariance = "i3")
# extract coefficients
est_coef_noNAs <- mirt.wrapper.coef(config_ord_noNAs)$coef
# Remove coefficients for i3 for HSLS
est_coef_noNAs[6, c("a1", "d1", "d2", "d3")] <- NA
# est_coef has columns group, item, a1, d1, d2, and d3. 
# a1 contains the loadings, d1-d3 are the thresholds. 
(els_coef_noNAs <- est_coef_noNAs[est_coef_noNAs$group == "ELS",])
(hsls_coef_noNAs <- est_coef_noNAs[est_coef_noNAs$group == "HSLS",])

# prepare lambda and nu matrices with estimates from the configural model
(lambda1_noNAs <- rbind(rep(els_coef_noNAs$a1, 3), rep(hsls_coef_noNAs$a1, 3)))
(nu1_noNAs <- rbind(unlist(els_coef_noNAs[, c("d1", "d2", "d3")]), 
              unlist(hsls_coef_noNAs[, c("d1", "d2", "d3")])))

# perform alignment and obtain aligned parameters of the latent mean and latent 
# variance
ord_align_noNAs <- invariance.alignment(lambda1_noNAs, nu1_noNAs, 
                                        wgt = sqrt(cbind(wgt_mat_noNAs, 
                                                         wgt_mat_noNAs,
                                                         wgt_mat_noNAs)))
ord_align_noNAs$pars
ord_align_noNAs$pars[2, 2]^2 # square the SD to get the variance, 0.955104

# Constrain latent means and variances to the aligned levels
mirtmodel_al_noNAs <- "
    F1 = i1, i2, i3, i4, i5
    START [HSLS] = (GROUP, MEAN_1, 0.6457948), (GROUP, COV_11, 0.955104)
  "
# Specify the mirt model explicitly providing the item names
mirtmodel_al_noNAs <- mirt.model(mirtmodel_al_noNAs, itemnames = m_items)

# Fit the aligned model
mod_aligned_noNAs <- multipleGroup(d_noNAs[, m_items], 
                             model = mirtmodel_al_noNAs, 
                             group = d_noNAs$sample, 
                             itemtype = "graded", 
                             invariance = "i3")
est_align_ord_noNAs <- mirt.wrapper.coef(mod_aligned_noNAs)$coef
est_align_ord_noNAs[6, c("a1", "d1", "d2", "d3")] <- NA

saveRDS(est_align_ord_noNAs, "rds/est_align_ord_noNAs.rds")
#  Compute factor scores and obtain standard errors of factor scores (ordinal case).
fs_align_ord_noNAs <- as.data.frame(fscores(mod_aligned_noNAs, mean = c(0, 0), cov = c(1, 1),
                                      full.scores.SE = TRUE, method = "EAP"))

# store the FS
score_df_noNAs <- as.data.frame(cbind(score_df_noNAs, 
    "partial_ord" = unlist(fs_partial_ord),    
    "approx_ord" = fs_align_ord_noNAs[, "F1"],
    "approx_ord_SE" = fs_align_ord_noNAs[, "SE_F1"]))

# Note that SE are not available for non-normal data in lavaan so we 
# cannot obtain SE for the partial invariance model assuming ordinal data.

score_df_noNAs[3:14] <- apply(score_df_noNAs[3:14], FUN = as.numeric, MARGIN = 2)
saveRDS(score_df_noNAs, "rds/score_df_noNAs.rds")



