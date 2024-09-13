
############# STEP 2: Traditional Measurement Invariance Testing ##############

# See `determine_configural_model.R` for details on choosing a configural model.

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
fit_config  <- cfa(cfa_config,  data = dat, group = "sample", estimator = "MLR",
                   missing = "FIML", se = "robust.mlr")
s_config <- summary(fit_config, fit.measures = TRUE)

###############################################################################
## STEP 2.1: Metric Invariance
###############################################################################

# Constrain loadings across groups, freely estimate intercepts and unique
# factor variances.

# Identification:
# Set the latent variance to 1 for the first group, free the latent variance in
# group 2.
# Set the latent means in both groups to 0.

# 2.1.1. Define a metric invariance model where all loadings are constrained 
# across groups
cfa_metric <-  'group: ELS
                math =~ l1 * NA * i1 + l2 * i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
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
                math =~ l1 * NA * i1 + l2 * i2 + l4 * i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric <- cfa(cfa_metric,  data = dat, estimator = "MLR",
                  group = "sample", missing = "FIML", se = "robust.mlr")
s_metric <- summary(fit_metric, fit.measures = TRUE)
LRT_conf_met <- lavTestLRT(fit_config, fit_metric)

# 2.1.2. Define a metric invariance model where one loading is freed across groups.

# 2.1.2.a. Free loading of excellentTests (i1)
metric1 <-     'group: ELS
                math =~ NA * i1 + l2 * i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
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
                math =~ NA * i1 + l2 * i2 + l4 * i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric1 <- cfa(metric1,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr")
s_metric1 <- summary(fit_metric1, fit.measures = TRUE)
LRT_met1 <- lavTestLRT(fit_metric, fit_metric1)

# 2.1.2.b. Free loading of understandTexts (i2)
metric2 <-     'group: ELS
                math =~ l1 * NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
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
                math =~ l1 * NA * i1 + i2 + l4 * i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric2 <- cfa(metric2,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr")
s_metric2 <- summary(fit_metric2, fit.measures = TRUE)
LRT_met2 <- lavTestLRT(fit_metric, fit_metric2)

# 2.1.2.c. Free loading of excellentAssign (i4)
metric4 <-     'group: ELS
                math =~ l1 * NA * i1 + l2 * i2 + l3 * i3 + i4 + l5 * i5
                       
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
                math =~ l1 * NA * i1 + l2 * i2 + i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric4 <- cfa(metric4,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr")
s_metric4 <- summary(fit_metric4, fit.measures = TRUE)
LRT_met4 <- lavTestLRT(fit_metric, fit_metric4)

# 2.1.2.d. Free loading of masterSkills (i5)
metric5 <-     'group: ELS
                math =~ l1 * NA * i1 + l2 * i2 + l3 * i3 + l4 * i4 + i5
                       
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
                math =~ l1 * NA * i1 + l2 * i2 + l4 * i4 + i5
                       
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
                math ~ 0 * 1  
                '
fit_metric5 <- cfa(metric5,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr",
                   std.lv = TRUE)
s_metric5 <- summary(fit_metric5, fit.measures = TRUE)
LRT_met5 <- lavTestLRT(fit_metric, fit_metric5)

# examine parameter estimates for the metric models considered
cbind(parameterEstimates(fit_metric)[,c(1:3)],
      "metric" = round(parameterEstimates(fit_metric)[,c(7)],3),
      "metric1" = round(parameterEstimates(fit_metric1)[,c(7)],3),
      "metric2" = round(parameterEstimates(fit_metric2)[,c(7)],3),
      "metric4" = round(parameterEstimates(fit_metric4)[,c(7)],3),
      "metric5" = round(parameterEstimates(fit_metric5)[,c(7)],3))

# fit statistics & LRT p values comparing the metric model to partial metric 
# models with one loading freed
options(scipen = 999) # turn off scientific notation (=0 for on)
fit_ind <- c("npar", "chisq", "df", "cfi", "tli", "rmsea", "srmr")
p_lrt_free1lambda <- c(LRT_conf_met$`Pr(>Chisq)`[2], 0, 
                       LRT_met1$`Pr(>Chisq)`[2], LRT_met2$`Pr(>Chisq)`[2], 
                       LRT_met4$`Pr(>Chisq)`[2], LRT_met5$`Pr(>Chisq)`[2])
p_lrt_free1chisqdiff <- c(LRT_conf_met$`Chisq diff`[2], 0,
                          LRT_met1$`Chisq diff`[2], LRT_met2$`Chisq diff`[2], 
                          LRT_met4$`Chisq diff`[2], LRT_met5$`Chisq diff`[2])
fit_tab_metric_1 <- round(rbind(cbind(s_config$fit[fit_ind],
                                      s_metric$fit[fit_ind], 
                                      s_metric1$fit[fit_ind],
                                      s_metric2$fit[fit_ind],
                                      s_metric4$fit[fit_ind],
                                      s_metric5$fit[fit_ind]),
                                p_lrt_free1lambda,
                                p_lrt_free1chisqdiff), 3)
colnames(fit_tab_metric_1) <- c("config", paste0("metric", c("", 1:2, 4:5)))
rownames(fit_tab_metric_1)[8] <- "p LRT"
rownames(fit_tab_metric_1)[9] <- "Chisq diff"
fit_tab_metric_1[8, 2] <- NA

saveRDS(fit_tab_metric_1, "rds/fit_tab_metric_1.rds")

# Freeing the loading of understandTexts (metric_2) results in the 
# largest improvement (largest chisq difference). We free the loading for 
# this item. 

# 2.1.3. Define a metric invariance model where two loadings are freed (loadings
# of i2 and one additional item).

# 2.1.3.a. Free loadings of understandTexts (i1) and excellentTests (i2)
metric_21 <-   'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
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
                math =~ NA * i1 + i2 + l4 * i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric_21 <- cfa(metric_21,  data = dat, estimator = "MLR",
                     group = "sample", missing = "FIML", se = "robust.mlr")
s_metric_21 <- summary(fit_metric_21, fit.measures = TRUE)
LRT_met21 <- lavTestLRT(fit_metric2, fit_metric_21)

saveRDS(fit_metric_21, "rds/fit_metric_21.rds")

# 2.1.3.b. Free loadings of excellentAssign (i2) and understandTexts (i4)
metric_24 <-   'group: ELS
                math =~ l1 * NA * i1 + i2 + l3 * i3 + i4 + l5 * i5
                       
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
                math =~ l1 * NA * i1 + i2 + i4 + l5 * i5
                       
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
                math ~ 0 * 1  
                '
fit_metric_24 <- cfa(metric_24,  data = dat, estimator = "MLR",
                     group = "sample", missing = "FIML", se = "robust.mlr",
                     std.lv = TRUE)
s_metric_24 <- summary(fit_metric_24, fit.measures = TRUE)
LRT_met24 <- lavTestLRT(fit_metric2, fit_metric_24)

# 2.1.3.c. Free loadings of understandTexts (i2) and masterSkills (i5)
metric_25 <-   'group: ELS
                math =~ l1 * NA * i1 + i2 + l3 * i3 + l4 * i4 + i5
                       
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
                math =~ l1 * NA * i1 + i2 + l4 * i4 + i5
                       
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
                math ~ 0 * 1  
                '
fit_metric_25 <- cfa(metric_25,  data = dat, estimator = "MLR",
                     group = "sample", missing = "FIML", se = "robust.mlr",
                     std.lv = TRUE)
s_metric_25 <- summary(fit_metric_25, fit.measures = TRUE)
LRT_met25 <- lavTestLRT(fit_metric2, fit_metric_25)

cbind(parameterEstimates(fit_metric)[,c(1:3)],
      "metric2" = round(parameterEstimates(fit_metric2)[,c(7)],3),
      "metric21" = round(parameterEstimates(fit_metric_21)[,c(7)],3),
      "metric24" = round(parameterEstimates(fit_metric_24)[,c(7)],3),
      "metric25" = round(parameterEstimates(fit_metric_25)[,c(7)],3))

# fit statistics & LRT p values comparing the metric model to partial metric 
# models with two loadings freed
p_lrt_free2lambda <- c(0, LRT_met21$`Pr(>Chisq)`[2],
                       LRT_met24$`Pr(>Chisq)`[2], LRT_met25$`Pr(>Chisq)`[2])
p_lrt_free1chisqdiff <- c(0, LRT_met21$`Chisq diff`[2], 
                          LRT_met24$`Chisq diff`[2],
                          LRT_met25$`Chisq diff`[2])
fit_tab_metric_2 <- round(rbind(cbind(s_metric2$fit[fit_ind],
                                      s_metric_21$fit[fit_ind],
                                      s_metric_24$fit[fit_ind],
                                      s_metric_25$fit[fit_ind]),
                                p_lrt_free2lambda, p_lrt_free1chisqdiff), 3)
colnames(fit_tab_metric_2) <- paste0("metric_2", c("", 1, 4:5))
rownames(fit_tab_metric_2)[8] <- "p LRT"
rownames(fit_tab_metric_2)[9] <- "Chisq diff"
fit_tab_metric_2[8:9,1] <- NA
#kableExtra::kable(fit_tab_metric_2)

saveRDS(fit_tab_metric_2, "rds/fit_tab_metric_2.rds")

# Freeing loadings for `excellentTests` and `understandTexts` (i2 and i1, in 
# metric_21) leads to the largest chi square difference. 

# We move on to examining partial scalar invariance with this model, as freeing 
# a third loading would lead the latent variable to be just identified in the 
# second group.

###############################################################################
## STEP 2.2: Scalar Invariance ## 
###############################################################################

# Additionally constrain the intercepts to be equal across groups. 

# Note that intercepts for `excellentTests` and `understandTexts` are freed along 
# with the loadings in line with Putnick & Bornstein (2016). This leaves two items 
# with invariant intercepts, `excellentAssign` and `masterSkills` (beyond 
# `understandClass`, which is missing in HSLS).
# Freeing either intercept would lead to a just-identified model, therefore we 
# do not test models with one intercept freed and only examine the scalar model.

# Identification:
# Setting the latent variance to 1 and the latent mean to 0 in the first group 
# is sufficient for identifiability. The latent mean and variance in the second 
# group are freely estimated.

# 2.2.1. Define a scalar invariance model where all intercepts are constrained 
# to be equal across groups
cfa_scalar <-  'group: ELS
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
fit_scalar <- cfa(cfa_scalar,  data = dat, estimator = "MLR",
                  group = "sample", missing = "FIML", se = "robust.mlr",
                  std.lv = TRUE)
s_scalar <- summary(fit_scalar, fit.measures = TRUE)
LRT_scalar <- lavTestLRT(fit_metric_21, fit_scalar)

saveRDS(fit_scalar, "rds/fit_scalar.rds")

###############################################################################
## STEP 2.3: Strict Invariance ## 
###############################################################################

# In addition to loadings and intercepts determined previously, constrain 
# all unique factor variances across groups.

# 2.3.1. Define a strict invariance model where all variances are constrained 
# to be equal across groups
cfa_strict <-  'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
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
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict <- cfa(cfa_strict,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr")
s_strict <- summary(fit_strict, fit.measures = TRUE)
LRT_strict <- lavTestLRT(fit_scalar, fit_strict)

# 2.3.2. Define a strict invariance model where one variance is freed.
# 2.3.2.a. Free variance of excellentTests (i1)
cfa_strict1 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
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
                i2 ~~ theta2 * i2
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict1 <- cfa(cfa_strict1,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr",
                   std.lv = TRUE)
s_strict1 <- summary(fit_strict1, fit.measures = TRUE)
LRT_strict1 <- lavTestLRT(fit_strict, fit_strict1)

# 2.3.2.b. Free variance of understandTexts (i2)
cfa_strict2 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2_1 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
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
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2_2 * i2
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict2 <- cfa(cfa_strict2,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr",
                   std.lv = TRUE)
s_strict2 <- summary(fit_strict2, fit.measures = TRUE)
LRT_strict2 <- lavTestLRT(fit_strict, fit_strict2)

# 2.3.2.c. Free variance of excellentAssign (i4)
cfa_strict4 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4_1 * i4
                i5 ~~ theta5 * i5
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
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i4 ~~ theta4_2 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict4 <- cfa(cfa_strict4,  data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr",
                   std.lv = TRUE)
s_strict4 <- summary(fit_strict4, fit.measures = TRUE)
(LRT_strict4 <- lavTestLRT(fit_strict, fit_strict4))

# 2.3.2.d. Free variance of masterSkills (i5) 
cfa_strict5 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
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
                       
                i1 ~~ theta1 * i1
                i2 ~~ theta2 * i2
                i4 ~~ theta4 * i4
                i5 ~~ theta5_2 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict5 <- cfa(cfa_strict4, data = dat, estimator = "MLR",
                   group = "sample", missing = "FIML", se = "robust.mlr",
                   std.lv = TRUE)
s_strict5 <- summary(fit_strict5, fit.measures = TRUE)
LRT_strict5 <- lavTestLRT(fit_strict, fit_strict5)

# fit statistics & LRT p values comparing the strict model to partial strict
# models with one variance freed
p_lrt_free1theta <- c(0, LRT_strict1$`Pr(>Chisq)`[2], 
                      LRT_strict2$`Pr(>Chisq)`[2],
                      LRT_strict4$`Pr(>Chisq)`[2], 
                      LRT_strict5$`Pr(>Chisq)`[2])
p_lrt_free1chisqdiff <- c(0, LRT_strict1$`Chisq diff`[2],
                          LRT_strict2$`Chisq diff`[2], 
                          LRT_strict4$`Chisq diff`[2],
                          LRT_strict5$`Chisq diff`[2])
fit_tab_strict <- round(rbind(cbind(s_strict$fit[fit_ind], 
                                    s_strict1$fit[fit_ind],
                                    s_strict2$fit[fit_ind], 
                                    s_strict4$fit[fit_ind],
                                    s_strict5$fit[fit_ind]),
                              p_lrt_free1theta, 
                              p_lrt_free1chisqdiff), 3)
colnames(fit_tab_strict) <- paste0("strict", c("", "1", "2", "4", "5"))
rownames(fit_tab_strict)[8] <- "p LRT"
rownames(fit_tab_strict)[9] <- "Chisq diff"
fit_tab_strict[8:9,1] <- NA
# Freeing the variance for i1 leads to the largest chisq difference.

saveRDS(fit_tab_strict, "rds/fit_tab_strict.rds")

# 2.3.3. Define a strict invariance model where two variances are freed,
# including that of i1.

# 2.3.3.a. Free variances of i2, i1
cfa_strict12 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2_1 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
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
                i4 ~~ theta4 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict12 <- cfa(cfa_strict12,  data = dat, estimator = "MLR",
                    group = "sample", missing = "FIML", se = "robust.mlr",
                    std.lv = TRUE)
s_strict12 <- summary(fit_strict12, fit.measures = TRUE)
LRT_strict12 <- lavTestLRT(fit_strict1, fit_strict12)

# 2.3.3.b. Free variances of i1, i4
cfa_strict14 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4_1 * i4
                i5 ~~ theta5 * i5
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
                i2 ~~ theta2 * i2
                i4 ~~ theta4_2 * i4
                i5 ~~ theta5 * i5
                i1 ~~ i2
                i2 ~~ i4
              
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict14 <- cfa(cfa_strict14,  data = dat, estimator = "MLR",
                    group = "sample", missing = "FIML", se = "robust.mlr",
                    std.lv = TRUE)
s_strict14 <- summary(fit_strict14, fit.measures = TRUE)
LRT_strict14 <- lavTestLRT(fit_strict1, fit_strict14)

# 2.3.3.a. Free variances of i1, i5
cfa_strict15 <- 'group: ELS
                math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                       
                i1 ~ NA * 1
                i2 ~ NA * 1
                i3 ~ nu3 * 1
                i4 ~ nu4 * 1
                i5 ~ nu5 * 1
                       
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4 * i4
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
                i2 ~~ theta2 * i2
                i4 ~~ theta4 * i4
                i5 ~~ theta5_2 * i5
                i1 ~~ i2
                i2 ~~ i4
                       
                math ~~ NA * math
                math ~ NA * 1  
                '
fit_strict15 <- cfa(cfa_strict15,  data = dat, estimator = "MLR",
                    group = "sample", missing = "FIML", se = "robust.mlr",
                    std.lv = TRUE)
s_strict15 <- summary(fit_strict15, fit.measures = TRUE)
LRT_strict15 <- lavTestLRT(fit_strict1, fit_strict15)

# fit statistics & LRT p values comparing the partial strict model with one
# variance freed to partial strict models with two variances freed
p_lrt_free2theta <- c(0, LRT_strict12$`Pr(>Chisq)`[2], 
                      LRT_strict14$`Pr(>Chisq)`[2],
                      LRT_strict15$`Pr(>Chisq)`[2])
p_lrt_free2chisqdiff <- c(0, LRT_strict12$`Chisq diff`[2], 
                          LRT_strict14$`Chisq diff`[2],
                          LRT_strict15$`Chisq diff`[2])
fit_tab_strict2 <- round(rbind(cbind(s_strict1$fit[fit_ind], 
                                     s_strict12$fit[fit_ind], 
                                     s_strict14$fit[fit_ind], 
                                     s_strict15$fit[fit_ind]),
                               p_lrt_free2theta, p_lrt_free2chisqdiff), 3)
colnames(fit_tab_strict2) <-  paste0("strict_1", c("", 2, 4:5))
rownames(fit_tab_strict2)[8] <- "p LRT"
rownames(fit_tab_strict2)[9] <- "Chisq diff"
fit_tab_strict2[8:9,1] <- NA

# Freeing the variance of item 2, understandTexts leads to the largest change in
# chisq difference.

saveRDS(fit_tab_strict2, "rds/fit_tab_strict2.rds")

# 2.3.4. Define a strict invariance model where three variances are freed, 
# including those for i1, i2.

# 2.3.4.a. Free variances of i2, i1, i4.
cfa_strict124 <- 'group: ELS
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
                  i5 ~~ theta5 * i5
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
                  i5 ~~ theta5 * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  '
fit_strict124 <- cfa(cfa_strict124,  data = dat, estimator = "MLR",
                     group = "sample", missing = "FIML", se = "robust.mlr",
                     std.lv = TRUE)
s_strict124 <- summary(fit_strict124, fit.measures = TRUE)
LRT_strict124 <- lavTestLRT(fit_strict12, fit_strict124)

saveRDS(fit_strict124, "rds/fit_strict124.rds")

# 2.3.4.a. Free variances of i2, i1, i5
cfa_strict125 <- 'group: ELS
                  math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                         
                  i1 ~ NA * 1
                  i2 ~ NA * 1
                  i3 ~ nu3 * 1
                  i4 ~ nu4 * 1
                  i5 ~ nu5 * 1
                         
                  i1 ~~ theta1_1 * i1
                  i2 ~~ theta2_1 * i2
                  i3 ~~ theta3 * i3
                  i4 ~~ theta4 * i4
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
                  i4 ~~ theta4 * i4
                  i5 ~~ theta5_2 * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  '
fit_strict125 <- cfa(cfa_strict125,  data = dat, estimator = "MLR",
                     group = "sample", missing = "FIML", se = "robust.mlr",
                     std.lv = TRUE)
s_strict125 <- summary(fit_strict125, fit.measures = TRUE)
LRT_strict125 <- lavTestLRT(fit_strict12, fit_strict125)

saveRDS(fit_strict125, "rds/fit_strict125.rds")

# fit statistics & LRT p values comparing the partial strict models with two vs. 
# three variances freed 
p_lrt_free3theta <- c(0, LRT_strict124$`Pr(>Chisq)`[2],
                      LRT_strict125$`Pr(>Chisq)`[2])
p_lrt_free3chisqdiff <- c(0, LRT_strict124$`Chisq diff`[2],
                          LRT_strict125$`Chisq diff`[2])
fit_tab_strict3 <- round(rbind(cbind(s_strict12$fit[fit_ind],
                                     s_strict124$fit[fit_ind],
                                     s_strict125$fit[fit_ind]),
                               p_lrt_free3theta, 
                               p_lrt_free3chisqdiff), 3)
colnames(fit_tab_strict3) <- c("strict_12", "strict_124","strict_125")
rownames(fit_tab_strict3)[8] <- "p LRT"
rownames(fit_tab_strict3)[9] <- "Chisq diff"
fit_tab_strict3[8:9,1] <- NA
# freeing the variance of i4 leads to the largest chisq difference.

saveRDS(fit_tab_strict3, "rds/fit_tab_strict3.rds")

# 2.3.5. Define a strict invariance model where four variances are freed.
# 2.3.5.a. Free variances of i2, i1, i4, i5
cfa_strict1245 <- 'group: ELS
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
fit_strict1245 <- cfa(cfa_strict1245,  data = dat, estimator = "MLR",
                      group = "sample", missing = "FIML", se = "robust.mlr",
                      std.lv = TRUE)
s_strict1245 <- summary(fit_strict1245, fit.measures = TRUE)
(LRT_strict1245 <- lavTestLRT(fit_strict124, fit_strict1245))

# cfa_strict1245 is our final partial invariance model.

saveRDS(fit_strict1245, "rds/fit_strict1245.rds")

# clean up environment
#rm(list=setdiff(ls(), c("dat", "fit_config", "s_config", "cfa_config")))
