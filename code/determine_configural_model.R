
############## STEP 1: Determine a Configural Model ###########################


### 1.0. ###
# define a basic model
cfa1 <- 'group: ELS
         math =~ i1 + i2 + i3 + i4 + i5
         
         group: HSLS
         math =~ i1 + i2 + i4 + i5
        '
# specify missing='FIML' to use full information maximum likelihood
cfa1_fit <- cfa(model = cfa1, data = dat, group = "sample", missing = "FIML", 
               estimator = "MLR", se = "robust.mlr")
summary(cfa1_fit, fit.measures = TRUE)

### 1.1. ###
# identify the model by freeing the first loading and setting the factor 
# variances to 1 and factor means to 0 in each group
cfa2 <- 'group: ELS
         math =~ NA * i1 + i2 + i3 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1
         
         group: HSLS
         math =~ NA * i1 + i2 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1        
        '
cfa2_fit <- cfa(model = cfa2, data = dat, group = "sample", missing = "FIML", 
                estimator = "MLR", se = "robust.mlr")
# can examine outputs from summary, lavCor
# summary(cfa2_fit, fit.measures = TRUE)
# lavCor(cfa2_fit)

### 1.2. ###
# examine modification indices
modindices(cfa2_fit, sort = TRUE)[1:10, c("lhs", "op", "rhs", "mi")]
# i1, i2 have the largest MI

### 1.3. ###
# add a correlation between i1, i2 in cfa2
cfa3 <- 'group: ELS
         math =~ NA * i1 + i2 + i3 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1
         i1 ~~ i2
         
         group: HSLS
         math =~ NA * i1 + i2 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1    
         i1 ~~ i2
        '
cfa3_fit <- cfa(model = cfa3, data = dat, group = "sample", missing = "FIML", 
                estimator = "MLR", se = "robust.mlr")
modindices(cfa3_fit, sort = TRUE)[1:10, c("lhs", "op", "rhs", "mi")]
# i2, i3 have the largest MI

### 1.4. ###
# update cfa3 with a correlation between i2, i3
cfa4 <- 'group: ELS
         math =~ NA * i1 + i2 + i3 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1
         i1 ~~ i2
         i2 ~~ i3
         
         group: HSLS
         math =~ NA * i1 + i2 + i4 + i5
         math ~~ 1 * math
         math ~ 0 * 1    
         i1 ~~ i2
        '
cfa4_fit <- cfa(model = cfa4, data = dat, group = "sample", missing = "FIML", 
                estimator = "MLR", se = "robust.mlr")
modindices(cfa4_fit, sort = TRUE)[1:10, c("lhs", "op", "rhs", "mi")]
# i1-i5, i2-i5, i2-i4 and i1-i4 have the largest MI. pick i2-i4 to correlate

### 1.5. ###
# update cfa4 with a correlation between i2, i4
cfa5 <-  'group: ELS
          math =~ NA * i1 + i2 + i3 + i4 + i5
          math ~~ 1 * math
          math ~ 0 * 1
          i1 ~~ i2
          i2 ~~ i3
          i2 ~~ i4
         
          group: HSLS
          math =~ NA * i1 + i2 + i4 + i5
          math ~~ 1 * math
          math ~ 0 * 1    
          i1 ~~ i2
          i2 ~~ i4
         '
cfa5_fit <- cfa(model = cfa5, data = dat, group = "sample", missing = "FIML", 
                estimator = "MLR", se = "robust.mlr")
modindices(cfa5_fit, sort = TRUE)[1:5, c("lhs", "op", "rhs", "mi")]
# we determine cfa5 as the final configural model


### final configural model ###
# redefine the configural model (cfa5) explicitly with more detailed labels:
cfa_config <-  'group: ELS
                math =~ NA * i1 + l2_1 * i2 + l3 * i3 + l4_1 * i4 + l5_1 * i5
                       
                i1 ~ nu1_1 * 1
                i2 ~ nu2_1 * 1
                i3 ~ nu3 * 1
                i4 ~ nu4_1 * 1
                i5 ~ nu5_1 * 1
                       
                # variances
                i1 ~~ theta1_1 * i1
                i2 ~~ theta2_1 * i2
                i3 ~~ theta3 * i3
                i4 ~~ theta4_1 * i4
                i5 ~~ theta5_1 * i5
                       
                # covariances
                i1 ~~ i2
                i2 ~~ cov3 * i3
                i2 ~~ i4
                       
                # latent variances
                math ~~ 1 * math
                # latent means
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

# clean up environment
rm(cfa1, cfa2, cfa2_fit, cfa3, cfa3_fit, cfa4, cfa4_fit, cfa5, cfa5_fit)
