############# STEP 4: Traditional Measurement Invariance Testing (ORDINAL) ##############

# Following tutorial by Tse, Lai, and Zhang: 
# Tse, W. W. Y., Lai, M. H., & Zhang, Y. (2024). Does strict invariance matter? 
# Valid group mean comparisons with ordered-categorical items. Behavior
# Research Methods, 56(4), 3117-3139.

###############################################################################
## STEP 4.0: Configural Invariance
###############################################################################

# We fit a single group one factor model to identify a marker variable.
cfa_m <-  'math =~ NA * i1 + i2 + i3 + i4 + i5
           math ~ 0 * 1 
           math ~~ 1 * math
'
parameterestimates(cfa(cfa_m, data = dat, ordered = TRUE, 
                       estimator = "WLSMV"))[1:7,]

# We choose item 4 as the marker as it has the largest loading estimate.

### 4.0.1. ###

# Constraints:
# Set the loading of the marker variable to be equal across groups (l4 in both).
# Set all intercepts to 0 in both groups.

# Identification:
# Fix the common factor variance to 1 for the first group and freely estimate 
# all loadings.
# Fix the common factor mean to 0 and the unique factor variances to 1 in the 
# first group.
# Fix one threshold for each item across groups, and two for the marker variable.
cfa_config_ord <- 'group: ELS
                   math =~ NA * i1 + i2 + i3 + l4 * i4 + i5
                     
                   # variances (set to 1 in the first group)
                   i1 ~~ 1 * i1
                   i2 ~~ 1 * i2
                   i3 ~~ 1 * i3
                   i4 ~~ 1 * i4
                   i5 ~~ 1 * i5
                  
                   # latent variance
                   math ~~ 1 * math
                   # latent mean
                   math ~ 0 * 1      
                  
                   # thresholds   
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
                   math =~ NA * i1 + i2 + l4 * i4 + i5
  
                   # unique factor variance (freely estimated in group 2)       
                   i1 ~~ NA * i1
                   i2 ~~ NA * i2
                   i4 ~~ NA * i4
                   i5 ~~ NA * i5
                         
                   # latent variance (freely est. in group 2)
                   math ~~ NA * math
                   # latent mean (freely est. in group 2)
                   math ~ NA * 1  
                  
                   # thresholds (initially set the first threshold for each item 
                   # equal across groups)
                   i1 | th1 * t1    
                   i1 | th2a * t2
                   i1 | th3a * t3
                  
                   i2 | th4 * t1    
                   i2 | th5a * t2
                   i2 | th6a * t3
                  
                   # for the marker variable, set an additional threshold equal
                   # across groups
                   i4 | th10 * t1    
                   i4 | th11 * t2
                   i4 | th12a * t3 
                  
                   i5 | th13 * t1    
                   i5 | th14a * t2
                   i5 | th15a * t3 
                  '
fit_config_ord  <- cfa(cfa_config_ord, data = dat, group = "sample", 
                       estimator = "WLSMV", ordered = TRUE,
                       parameterization = "theta")
s_config_ord <- summary(fit_config_ord, fit.measures = TRUE)

# Examine modification indices
modindices(fit_config_ord, sort = TRUE)[1:5, c("lhs", "op", "rhs", "mi")]
# Add a covariance between i1 and i2 as this item has the largest modification index.

### 4.0.2. ###
cfa_config_ord_c <-  'group: ELS
                  math =~ NA * i1 + i2 + i3 + l4 * i4 + i5
                     
                  i1 ~~ 1 * i1
                  i2 ~~ 1 * i2
                  i3 ~~ 1 * i3
                  i4 ~~ 1 * i4
                  i5 ~~ 1 * i5
                  i1 ~~ i2
                  
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
                  math =~ NA * i1 + i2 + l4 * i4 + i5
  
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_config_ord_c  <- cfa(cfa_config_ord_c, data = dat, group = "sample", 
                         estimator = "WLSMV", ordered = TRUE,
                         parameterization = "theta")
s_config_ord_c <- summary(fit_config_ord_c, fit.measures = TRUE)

modindices(fit_config_ord_c, sort = TRUE)[1:5, c("lhs", "op", "rhs", "mi")]
# For ELS, we add a covariance between i2 and i3.

### 4.0.3. ###
cfa_config_ord_c2 <- 'group: ELS
                      math =~ NA * i1 + i2 + i3 + l4 * i4 + i5
                         
                      i1 ~~ 1 * i1
                      i2 ~~ 1 * i2
                      i3 ~~ 1 * i3
                      i4 ~~ 1 * i4
                      i5 ~~ 1 * i5
                      i1 ~~ i2
                      i2 ~~ i3
                      
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
                      math =~ NA * i1 + i2 + l4 * i4 + i5
      
                      i1 ~~ NA * i1
                      i2 ~~ NA * i2
                      i4 ~~ NA * i4
                      i5 ~~ NA * i5
                      i1 ~~ i2
                             
                      math ~~ NA * math
                      math ~ NA * 1  
                      
                      i1 | th1 * t1    
                      i1 | th2a * t2
                      i1 | th3a * t3
                      
                      i2 | th4 * t1    
                      i2 | th5a * t2
                      i2 | th6a * t3
                      
                      i4 | th10 * t1    
                      i4 | th11 * t2
                      i4 | th12a * t3 
                      
                      i5 | th13 * t1    
                      i5 | th14a * t2
                      i5 | th15a * t3 
                      '
fit_config_ord_c2  <- cfa(cfa_config_ord_c2, data = dat, group = "sample", 
                          estimator = "WLSMV", ordered = TRUE,
                          parameterization = "theta")
s_config_ord_c2 <- summary(fit_config_ord_c2, fit.measures = TRUE)

modindices(fit_config_ord_c2, sort = TRUE)[1:5, c("lhs", "op", "rhs", "mi")]

# i1-i4, i1-i5, i2-i5 and i2-i4 have the same and the largest modification index. 
# We choose to add a correlation between i2 and i4.

### 4.0.4. ### 
cfa_config_ord_c3 <- 'group: ELS
                      math =~ NA * i1 + i2 + i3 + l4 * i4 + i5
                         
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
                      math =~ NA * i1 + i2 + l4 * i4 + i5
      
                      i1 ~~ NA * i1
                      i2 ~~ NA * i2
                      i4 ~~ NA * i4
                      i5 ~~ NA * i5
                      i1 ~~ i2
                      i2 ~~ i4
                      
                      math ~~ NA * math
                      math ~ NA * 1  
                      
                      i1 | th1 * t1    
                      i1 | th2a * t2
                      i1 | th3a * t3
                      
                      i2 | th4 * t1    
                      i2 | th5a * t2
                      i2 | th6a * t3
                      
                      i4 | th10 * t1    
                      i4 | th11 * t2
                      i4 | th12a * t3 
                      
                      i5 | th13 * t1    
                      i5 | th14a * t2
                      i5 | th15a * t3 
                      '
fit_config_ord_c3  <- cfa(cfa_config_ord_c3, data = dat, group = "sample", 
                          estimator = "WLSMV", ordered = TRUE,
                          parameterization = "theta")
s_config_ord_c3 <- summary(fit_config_ord_c3, fit.measures = TRUE)

modindices(fit_config_ord_c3, sort = TRUE)[1:5, c("lhs", "op", "rhs", "mi")]

# We proceed with `fit_config_ord_c3` as our final configural model.


###############################################################################
## STEP 4.1: Metric Invariance
###############################################################################

### 4.1.0. Define a metric invariance model where all loadings are constrained 
# across groups
cfa_metric <-    'group: ELS
                  math =~ NA * l1 * i1 + l2 * i2 + l3 * i3 + l4 * i4 + l5 * i5
                     
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
                  
                  # thresholds   
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
                  math =~ NA * l1 * i1 + l2 * i2 + l4 * i4 + l5 * i5
    
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric  <- cfa(cfa_metric, data = dat, group = "sample", 
                   estimator = "WLSMV", ordered = TRUE,
                   parameterization = "theta")
s_metric <- summary(fit_metric, fit.measures = TRUE)
# The metric model has a scaled RMSEA of 0.049, which indicates acceptable fit 
# and SRMR of 0.003, which indicates good fit.

(LRT_conf4_met <- lavTestLRT(fit_config_ord_c3, fit_metric))
# The chi-square test indicates that the metric model has a significantly worse 
# fit. We proceed to determining a partial metric model by releasing one loading.

### 4.1.1. Define a metric invariance model where one loading is freed across groups.

# 4.1.1.a. Free loading of excellentTests (i1)
cfa_metric_1 <-  'group: ELS
                  math =~ NA * i1 + l2 * i2 + l3 * i3 + l4 * i4 + l5 * i5
                     
                  # variances (set to 1 in the first group)
                  i1 ~~ 1 * i1
                  i2 ~~ 1 * i2
                  i3 ~~ 1 * i3
                  i4 ~~ 1 * i4
                  i5 ~~ 1 * i5
                  i1 ~~ i2
                  i2 ~~ i3
                  i2 ~~ i4
                  
                  # latent variance
                  math ~~ 1 * math
                  # latent mean
                  math ~ 0 * 1      
                  
                  # thresholds   
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
                  math =~ NA * i1 + l2 * i2 + l4 * i4 + l5 * i5
  
                  # unique factor variance (freely estimated in group 2)       
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  # latent variance (freely est. in group 2)
                  math ~~ NA * math
                  # latent mean (freely est. in group 2)
                  math ~ NA * 1  
                  
                  # thresholds (initially set the first threshold for each item 
                  # equal across groups)
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  # for the marker variable, set an additional threshold equal
                  # across groups
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_1  <- cfa(cfa_metric_1, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_metric_1 <- summary(fit_metric_1, fit.measures = TRUE)

(LRT_met_met_1 <- lavTestLRT(fit_metric, fit_metric_1))

# 4.1.1.b. Free loading of understandTexts (i2)
cfa_metric_2 <-  'group: ELS
                  math =~ NA * l1 * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
                     
                  # variances (set to 1 in the first group)
                  i1 ~~ 1 * i1
                  i2 ~~ 1 * i2
                  i3 ~~ 1 * i3
                  i4 ~~ 1 * i4
                  i5 ~~ 1 * i5
                  i1 ~~ i2
                  i2 ~~ i3     
                  i2 ~~ i4
                  
                  # latent variance
                  math ~~ 1 * math
                  # latent mean
                  math ~ 0 * 1      
                  
                  # thresholds   
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
                  math =~ NA * l1 * i1 + i2 + l4 * i4 + l5 * i5
        
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_2  <- cfa(cfa_metric_2, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_metric_2 <- summary(fit_metric_2, fit.measures = TRUE)

(LRT_met_met_2 <- lavTestLRT(fit_metric, fit_metric_2))

# 4.1.1.c. Free loading of excellentAssign (i4)
cfa_metric_4 <-  'group: ELS
                  math =~ NA * l1 * i1 + l2 * i2 + l3 * i3 + i4 + l5 * i5
                     
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
                  math =~ NA * l1 * i1 + l2 * i2 + i4 + l5 * i5
  
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_4  <- cfa(cfa_metric_4, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_metric_4 <- summary(fit_metric_4, fit.measures = TRUE)

(LRT_met_met_4 <- lavTestLRT(fit_metric, fit_metric_4))

# 4.1.1.d. Free loading of masterSkills (i5)
cfa_metric_5 <-  'group: ELS
                  math =~ NA * l1 * i1 + l2 * i2 + l3 * i3 + l4 * i4 + i5
                     
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
                  math =~ NA * l1 * i1 + l2 * i2 + l4 * i4 + i5
  
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                         
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_5  <- cfa(cfa_metric_5, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_metric_5 <- summary(fit_metric_5, fit.measures = TRUE)

(LRT_met_met_5 <- lavTestLRT(fit_metric, fit_metric_5))


# fit statistics & LRT p values comparing the metric model to partial metric 
# models with one loading freed
options(scipen = 999) # turn off scientific notation (=0 for on)
fit_ind <- c("npar", "chisq", "df", "cfi", "tli", "rmsea.scaled", "srmr")
p_lrt_free1lambda <- c(LRT_conf4_met$`Pr(>Chisq)`[2], 0,
                       LRT_met_met_1$`Pr(>Chisq)`[2],
                       LRT_met_met_2$`Pr(>Chisq)`[2], 
                       LRT_met_met_4$`Pr(>Chisq)`[2],
                       LRT_met_met_5$`Pr(>Chisq)`[2])
p_lrt_free1chisqdiff <- c(LRT_conf4_met$`Chisq diff`[2], 0,
                          LRT_met_met_1$`Chisq diff`[2],
                          LRT_met_met_2$`Chisq diff`[2], 
                          LRT_met_met_4$`Chisq diff`[2],
                          LRT_met_met_5$`Chisq diff`[2])
fit_tab_metric_1 <- round(rbind(cbind(s_config_ord_c3$fit[fit_ind],
                                      s_metric$fit[fit_ind], 
                                      s_metric_1$fit[fit_ind],
                                      s_metric_2$fit[fit_ind],
                                      s_metric_4$fit[fit_ind],
                                      s_metric_5$fit[fit_ind]),
                                p_lrt_free1lambda,
                                p_lrt_free1chisqdiff), 3)
colnames(fit_tab_metric_1) <- c("config", paste0("metric", c("", 1:2, 4:5)))
rownames(fit_tab_metric_1)[8] <- "p LRT"
rownames(fit_tab_metric_1)[9] <- "Chisq diff"
fit_tab_metric_1[8, 2] <- NA
fit_tab_metric_1


# The chisq difference is largest for `metric2`, which suggests we should free 
# the loading for i2 next. As there is insufficient evidence that item 4, 
# which we had selected as our marker variable, is noninvariant, we may proceed
# with our analyses in line with Tse, Lai, and Zhang's (2024) suggestions.

### 4.1.2. Define a metric invariance model where two loadings are freed 
# (loadings of i2 and one additional item).
 
# 4.1.2.a. Free loadings of understandTexts (i2) and excellentTests (i1)
cfa_metric_21 <- 'group: ELS
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
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_21  <- cfa(cfa_metric_21, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_metric_21 <- summary(fit_metric_21, fit.measures = TRUE)

(LRT_met_2_met_21 <- lavTestLRT(fit_metric_2, fit_metric_21))


# 4.1.2.b. Free loadings of  excellentAssign (i4) and understandTexts (i2)
cfa_metric_24 <- 'group: ELS
                  math =~ NA * l1 * i1 + i2 + l3 * i3 + i4 + l5 * i5
                  
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
                  math =~ NA * l1 * i1 + i2 + i4 + l5 * i5
                  
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_24  <- cfa(cfa_metric_24, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_metric_24 <- summary(fit_metric_24, fit.measures = TRUE)

(LRT_met_2_met_24 <- lavTestLRT(fit_metric_2, fit_metric_24))

# 4.1.2.c. Free loadings of understandTexts (i2) and masterSkills (i5)
cfa_metric_25 <- 'group: ELS
                  math =~ NA * l1 * i1 + i2 + l3 * i3 + l4 * i4 + i5
                  
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
                  math =~ NA * l1 * i1 + i2 + l4 * i4 + i5
                  
                  i1 ~~ NA * i1
                  i2 ~~ NA * i2
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3a * t3
                  
                  i2 | th4 * t1    
                  i2 | th5a * t2
                  i2 | th6a * t3
                  
                  i4 | th10 * t1    
                  i4 | th11 * t2
                  i4 | th12a * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15a * t3 
                  '
fit_metric_25  <- cfa(cfa_metric_25, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_metric_25 <- summary(fit_metric_25, fit.measures = TRUE)

(LRT_met_2_met_25 <- lavTestLRT(fit_metric_2, fit_metric_25))

# fit statistics & LRT p values comparing the metric model to partial metric 
# models with one loading freed
options(scipen = 999) # turn off scientific notation (=0 for on)
fit_ind <- c("npar", "chisq", "df", "cfi", "tli", "rmsea.scaled", "srmr")
p_lrt_free2lambda <- c(LRT_met_met_2$`Pr(>Chisq)`[2], 0,
                        LRT_met_2_met_21$`Pr(>Chisq)`[2], 
                        LRT_met_2_met_24$`Pr(>Chisq)`[2],
                        LRT_met_2_met_25$`Pr(>Chisq)`[2])
p_lrt_free2chisqdiff <- c(LRT_met_met_2$`Chisq diff`[2], 0,
                           LRT_met_2_met_21$`Chisq diff`[2], 
                           LRT_met_2_met_24$`Chisq diff`[2],
                           LRT_met_2_met_25$`Chisq diff`[2])
fit_tab_metric_2 <- round(rbind(cbind(s_metric$fit[fit_ind],
                                       s_metric_2$fit[fit_ind],
                                       s_metric_21$fit[fit_ind],
                                       s_metric_24$fit[fit_ind],
                                       s_metric_25$fit[fit_ind]),
                                 p_lrt_free2lambda, p_lrt_free2chisqdiff), 3)
colnames(fit_tab_metric_2) <- c("metric", paste0("metric_2", c("", 1, 4:5)))
rownames(fit_tab_metric_2)[8] <- "p LRT"
rownames(fit_tab_metric_2)[9] <- "Chisq diff"
fit_tab_metric_2[8, 2] <- NA
fit_tab_metric_2


# The largest $\chi^2$ difference is observed for `metric_21`. We proceed with 
# the model where the loadings are freed for items 2 and 1 to test for scalar 
# invariance, as releasing the loading of i4 or i5 would lead to the latent 
# variable being just identified in the second group.

###############################################################################
## STEP 4.2: Scalar Invariance
###############################################################################

### 4.2.0. Define a scalar invariance model where all thresholds are constrained 
# across groups
cfa_scalar <-'group: ELS
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
              i4 ~~ NA * i4
              i5 ~~ NA * i5
              i1 ~~ i2
              i2 ~~ i4
              
              math ~~ NA * math
              math ~ NA * 1  
              
              i1 | th1 * t1    
              i1 | th2 * t2
              i1 | th3 * t3
              
              i2 | th4 * t1    
              i2 | th5 * t2
              i2 | th6 * t3
              
              i4 | th10 * t1    
              i4 | th11 * t2
              i4 | th12 * t3 
              
              i5 | th13 * t1    
              i5 | th14 * t2
              i5 | th15 * t3 
              '
fit_scalar  <- cfa(cfa_scalar, data = dat, group = "sample", 
                   estimator = "WLSMV", ordered = TRUE,
                   parameterization = "theta")
s_scalar <- summary(fit_scalar, fit.measures = TRUE)

(LRT_met_21_scalar <- lavTestLRT(fit_metric_21, fit_scalar))

modindices(fit_scalar, free.remove = "FALSE", op = "|", sort = TRUE)

### 4.1.0.a. Release the second threshold for i1.
cfa_s_i1t2 <-'group: ELS
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
              i4 ~~ NA * i4
              i5 ~~ NA * i5
              i1 ~~ i2
              i2 ~~ i4
              
              math ~~ NA * math
              math ~ NA * 1  
              
              i1 | th1 * t1    
              i1 | th2a * t2
              i1 | th3 * t3
              
              i2 | th4 * t1    
              i2 | th5 * t2
              i2 | th6 * t3
              
              i4 | th10 * t1    
              i4 | th11 * t2
              i4 | th12 * t3 
              
              i5 | th13 * t1    
              i5 | th14 * t2
              i5 | th15 * t3 
              '
fit_s_i1t2  <- cfa(cfa_s_i1t2, data = dat, group = "sample", 
                   estimator = "WLSMV", ordered = TRUE,
                   parameterization = "theta")
s_s_i1t2 <- summary(fit_s_i1t2, fit.measures = TRUE)

(LRT_s_i1t2 <- lavTestLRT(fit_scalar, fit_s_i1t2))

modindices(fit_s_i1t2, free.remove = "FALSE", op = "|", sort = TRUE)[1:10,]

### 4.1.0.b. Release the second threshold for i1 and the second threshold of i4.
cfa_s_i1t2_i4t2 <- 'group: ELS
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
                    i4 ~~ NA * i4
                    i5 ~~ NA * i5
                    i1 ~~ i2
                    i2 ~~ i4
                    
                    math ~~ NA * math
                    math ~ NA * 1  
                    
                    i1 | th1 * t1    
                    i1 | th2a * t2
                    i1 | th3 * t3
                    
                    i2 | th4 * t1    
                    i2 | th5 * t2
                    i2 | th6 * t3
                    
                    i4 | th10 * t1    
                    i4 | th11a * t2
                    i4 | th12 * t3 
                    
                    i5 | th13 * t1    
                    i5 | th14 * t2
                    i5 | th15 * t3 
                    '
fit_s_i1t2_i4t2  <- cfa(cfa_s_i1t2_i4t2, data = dat, group = "sample", 
                        estimator = "WLSMV", ordered = TRUE,
                        parameterization = "theta")
s_s_i1t2_i4t2 <- summary(fit_s_i1t2_i4t2, fit.measures = TRUE)

(LRT_s_i1t2_i4t2 <- lavTestLRT(fit_s_i1t2, fit_s_i1t2_i4t2))

modindices(fit_s_i1t2_i4t2, 
           free.remove = "FALSE", op = "|", sort = TRUE)[1:10,]

### 4.1.0.c. Release the second threshold for i1, the second threshold of i4, 
# and the second threshold of i5.
cfa_s_i1t2_i4t2_i5t2 <- 
                 'group: ELS
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
                  i4 ~~ NA * i4
                  i5 ~~ NA * i5
                  i1 ~~ i2
                  i2 ~~ i4
                  
                  math ~~ NA * math
                  math ~ NA * 1  
                  
                  i1 | th1 * t1    
                  i1 | th2a * t2
                  i1 | th3 * t3
                  
                  i2 | th4 * t1    
                  i2 | th5 * t2
                  i2 | th6 * t3
                  
                  i4 | th10 * t1    
                  i4 | th11a * t2
                  i4 | th12 * t3 
                  
                  i5 | th13 * t1    
                  i5 | th14a * t2
                  i5 | th15 * t3 
                  '
fit_s_i1t2_i4t2_i5t2  <- cfa(cfa_s_i1t2_i4t2_i5t2,  data = dat, 
                             group = "sample", estimator = "WLSMV", 
                             ordered = TRUE, parameterization = "theta")
s_s_i1t2_i4t2_i5t2 <- summary(fit_s_i1t2_i4t2_i5t2, fit.measures = TRUE)

(LRT_sca_i1t2_i4t2_i5t2 <- lavTestLRT(fit_s_i1t2_i4t2, fit_s_i1t2_i4t2_i5t2))

modindices(fit_s_i1t2_i4t2_i5t2, 
           free.remove = "FALSE", op = "|", sort = TRUE)[1:10,]

### 4.1.0.d. Release the second threshold for i1, the second threshold of i4, 
# the second threshold of i5, and the second threshold of i2.
cfa_s_i1t2_i4t2_i5t2_i2t2 <- 
                 'group: ELS
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
                  i4 ~~ NA * i4
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
                  i5 | th15 * t3 
                  '
fit_s_i1t2_i4t2_i5t2_i2t2  <- cfa(cfa_s_i1t2_i4t2_i5t2_i2t2,  data = dat,
                                  group = "sample", 
                                  estimator = "WLSMV", ordered = TRUE,
                                  parameterization = "theta")
s_s_i1t2_i4t2_i5t2_i2t2 <- summary(fit_s_i1t2_i4t2_i5t2_i2t2, fit.measures = TRUE)

(LRT_sca_i1t2_i4t2_i5t2_i2t2 <- 
   lavTestLRT(fit_s_i1t2_i4t2_i5t2, fit_s_i1t2_i4t2_i5t2_i2t2))

modindices(fit_s_i1t2_i4t2_i5t2_i2t2, 
           free.remove = "FALSE", op = "|", sort = TRUE)[1:10,]

### 4.1.0.e. Release the second threshold for i1, the second threshold of i4, 
# the second threshold of i5, the second threshold of i2, and the third 
# threshold of i5.
cfa_s_i1t2_i4t2_i5t2_i2t2_i5i3 <- 
                 'group: ELS
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
                  i4 ~~ NA * i4
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
fit_s_i1t2_i4t2_i5t2_i2t2_i5i3  <- cfa(cfa_s_i1t2_i4t2_i5t2_i2t2_i5i3,  
                                       data = dat, group = "sample", 
                                       estimator = "WLSMV", ordered = TRUE,
                                       parameterization = "theta")
s_s_i1t2_i4t2_i5t2_i2t2_i5i3 <- summary(fit_s_i1t2_i4t2_i5t2_i2t2_i5i3, 
                                   fit.measures = TRUE)

(LRT_sca_i1t2_i4t2_i5t2_i2t2_i5t3 <- 
   lavTestLRT(fit_s_i1t2_i4t2_i5t2_i2t2, fit_s_i1t2_i4t2_i5t2_i2t2_i5i3))

modindices(fit_s_i1t2_i4t2_i5t2_i2t2_i5i3, 
           free.remove = "FALSE", op = "|", sort = TRUE)[1:10,]
# As the modification indices have fallen below 3.84, we release no more 
# thresholds and proceed to testing strict invariance by constraining all 
# unique variances across groups.

###############################################################################
## STEP 4.3: Strict Invariance
###############################################################################

### 4.3.0. Define a strict invariance model where all variances to are set to 
# one in both groups.
cfa_strict <-  'group: ELS
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
                
                i1 ~~ 1 * i1
                i2 ~~ 1 * i2
                i4 ~~ 1 * i4
                i5 ~~ 1 * i5
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
fit_strict  <- cfa(cfa_strict, data = dat, group = "sample", 
                   estimator = "WLSMV", ordered = TRUE,
                   parameterization = "theta")
s_strict <- summary(fit_strict, fit.measures = TRUE)

(LRT_partial_scalar_strict <- 
   lavTestLRT(fit_s_i1t2_i4t2_i5t2_i2t2_i5i3, fit_strict))


# 4.3.1. Define a strict invariance model where one variance is freed.

# 4.3.1.a. Free variance of excellentTests (i1) in group 2
cfa_strict_1 <-  'group: ELS
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
                  i2 ~~ 1 * i2
                  i4 ~~ 1 * i4
                  i5 ~~ 1 * i5
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
fit_strict_1  <- cfa(cfa_strict_1, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_strict_1 <- summary(fit_strict_1, fit.measures = TRUE)


(LRT_strict_strict1 <- lavTestLRT(fit_strict, fit_strict_1))

# 4.3.1.b. Free variance of understandTexts (i2) in group 2
cfa_strict_2 <-  'group: ELS
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
                  
                  i1 ~~ 1 * i1
                  i2 ~~ NA * i2
                  i4 ~~ 1 * i4
                  i5 ~~ 1 * i5
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
fit_strict_2  <- cfa(cfa_strict_2, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_strict_2 <- summary(fit_strict_2, fit.measures = TRUE)

(LRT_strict_strict2 <- lavTestLRT(fit_strict, fit_strict_2))

# 4.3.1.c. Free variance of excellentAssign (i4) in group 2
cfa_strict_4 <-  'group: ELS
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
                  
                  i1 ~~ 1 * i1
                  i2 ~~ 1 * i2
                  i4 ~~ NA * i4
                  i5 ~~ 1 * i5
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
fit_strict_4  <- cfa(cfa_strict_4, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_strict_4 <- summary(fit_strict_4, fit.measures = TRUE)

(LRT_strict_strict4 <- lavTestLRT(fit_strict, fit_strict_4))

# 4.3.1.d. Free variance of masterSkills (i5) in group 2
cfa_strict_5 <-  'group: ELS
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
                  
                  i1 ~~ 1 * i1
                  i2 ~~ 1 * i2
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
fit_strict_5  <- cfa(cfa_strict_5, data = dat, group = "sample", 
                     estimator = "WLSMV", ordered = TRUE,
                     parameterization = "theta")
s_strict_5 <- summary(fit_strict_5, fit.measures = TRUE)

(LRT_strict_strict5 <- lavTestLRT(fit_strict, fit_strict_5))

# fit statistics & LRT p values comparing the metric model to partial metric 
# models with one loading freed
options(scipen = 999) # turn off scientific notation (=0 for on)
fit_ind <- c("npar", "chisq", "df", "cfi", "tli", "rmsea.scaled", "srmr")
p_lrt_free1theta <- c(LRT_partial_scalar_strict$`Pr(>Chisq)`[2], 0,
                       LRT_strict_strict1$`Pr(>Chisq)`[2],
                       LRT_strict_strict2$`Pr(>Chisq)`[2], 
                       LRT_strict_strict4$`Pr(>Chisq)`[2],
                       LRT_strict_strict5$`Pr(>Chisq)`[2])
p_lrt_free1chisqdiff_th <- c(LRT_partial_scalar_strict$`Chisq diff`[2], 0,
                          LRT_strict_strict1$`Chisq diff`[2],
                          LRT_strict_strict2$`Chisq diff`[2], 
                          LRT_strict_strict4$`Chisq diff`[2],
                          LRT_strict_strict5$`Chisq diff`[2])
fit_tab_strict_1 <- round(rbind(cbind(s_s_i1t2_i4t2_i5t2_i2t2_i5i3$fit[fit_ind],
                                      s_strict$fit[fit_ind], 
                                      s_strict_1$fit[fit_ind],
                                      s_strict_2$fit[fit_ind],
                                      s_strict_4$fit[fit_ind],
                                      s_strict_5$fit[fit_ind]),
                                p_lrt_free1theta,
                                p_lrt_free1chisqdiff_th), 3)
colnames(fit_tab_strict_1) <- c("partial scalar", paste0("strict", c("", 1:2, 4:5)))
rownames(fit_tab_strict_1)[8] <- "p LRT"
rownames(fit_tab_strict_1)[9] <- "Chisq diff"
fit_tab_strict_1[8, 2] <- NA
fit_tab_strict_1

# Release the unique variance for i1, proceed with testing models with two
# unique variances released.

# 4.3.2. Define a strict invariance model where two variances are freed.

# 4.3.2.a. Free variance of i1, i2
cfa_strict_12 <-   'group: ELS
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
                    i5 ~~ 1 * i5
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
fit_strict_12  <- cfa(cfa_strict_12, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_strict_12 <- summary(fit_strict_12, fit.measures = TRUE)

(LRT_strict1_strict12 <- lavTestLRT(fit_strict_1, fit_strict_12))

# 4.3.2.b. Free variances of i1, i4 in group 2
cfa_strict_14 <- 'group: ELS
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
                  i2 ~~ 1 * i2
                  i4 ~~ NA * i4
                  i5 ~~ 1 * i5
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
fit_strict_14  <- cfa(cfa_strict_14, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_strict_14 <- summary(fit_strict_14, fit.measures = TRUE)

(LRT_strict1_strict14 <- lavTestLRT(fit_strict_1, fit_strict_14))

# 4.3.2.c. Free variances of i1, i5 in group 2
cfa_strict_15 <-   'group: ELS
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
                    i2 ~~ 1 * i2
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
fit_strict_15  <- cfa(cfa_strict_15, data = dat, group = "sample", 
                      estimator = "WLSMV", ordered = TRUE,
                      parameterization = "theta")
s_strict_15 <- summary(fit_strict_15, fit.measures = TRUE)

(LRT_strict1_strict15 <- lavTestLRT(fit_strict_1, fit_strict_15))

# fit statistics & LRT p values comparing the metric model to partial metric 
# models with one loading freed
options(scipen = 999) # turn off scientific notation (=0 for on)
fit_ind <- c("npar", "chisq", "df", "cfi", "tli", "rmsea.scaled", "srmr")
p_lrt_free2theta <- c(LRT_strict_strict1$`Pr(>Chisq)`[2], 0,
                        LRT_strict1_strict12$`Pr(>Chisq)`[2], 
                        LRT_strict1_strict14$`Pr(>Chisq)`[2],
                        LRT_strict1_strict15$`Pr(>Chisq)`[2])
p_lrt_free2chisqdiff_th <- c(LRT_strict_strict1$`Chisq diff`[2], 0,
                           LRT_strict1_strict12$`Chisq diff`[2], 
                           LRT_strict1_strict14$`Chisq diff`[2],
                           LRT_strict1_strict15$`Chisq diff`[2])
fit_tab_strict_2 <- round(rbind(cbind(s_strict$fit[fit_ind],
                                       s_strict_1$fit[fit_ind],
                                       s_strict_12$fit[fit_ind],
                                       s_strict_14$fit[fit_ind],
                                       s_strict_15$fit[fit_ind]),
                                 p_lrt_free2theta, p_lrt_free2chisqdiff_th), 3)
colnames(fit_tab_strict_2) <- c("strict", paste0("strict_2", c("", 1, 4:5)))
rownames(fit_tab_strict_2)[8] <- "p LRT"
rownames(fit_tab_strict_2)[9] <- "Chisq diff"
fit_tab_strict_2[8, 2] <- NA
fit_tab_strict_2

# $\chi^2$ difference is largest for `strict_12`. We test versions of this model 
# with three variances released.

# 4.3.3. Define a strict invariance model where three variances are freed.

# 4.3.3.a. Free variances of i1, i2, i4 in group 2.
cfa_strict_124 <-'group: ELS
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
fit_strict_124  <- cfa(cfa_strict_124, data = dat, group = "sample", 
                       estimator = "WLSMV", ordered = TRUE,
                       parameterization = "theta")
s_strict_124 <- summary(fit_strict_124, fit.measures = TRUE)
#lavInspect(fit_strict_124, what = "est")
saveRDS(cfa_strict_124, "rds/cfa_partial_ord.rds")
saveRDS(fit_strict_124, "rds/fit_partial_ord.rds")

(LRT_strict12_strict124 <- lavTestLRT(fit_strict_12, fit_strict_124))

# 4.3.3.b. Free variances of i1, i2, i5 in group 2.
cfa_strict_125 <-  'group: ELS
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
fit_strict_125  <- cfa(cfa_strict_125, data = dat, group = "sample", 
                       estimator = "WLSMV", ordered = TRUE,
                       parameterization = "theta")
s_strict_125 <- summary(fit_strict_125, fit.measures = TRUE)

(LRT_strict12_strict125 <- lavTestLRT(fit_strict_12, fit_strict_125))
# Both LRTs are nonsignificant, with the same $\chi^2$ difference. We decide to 
# release the variance for i4. The final partial strict invariance model is 
# `fit_strict_124`.

# Factor scores

fit_partial_ord <- fit_strict_124

fit_partial_ord  <- cfa(cfa_strict_124, data = dat, group = "sample", 
                        estimator = "WLSMV", ordered = TRUE,
                        parameterization = "theta")

fs_partial_ord <- lavPredict(fit_partial_ord, method = "EBM", se = TRUE)

# saveRDS(fs_partial_ord, "fs_partial_ord.rds")