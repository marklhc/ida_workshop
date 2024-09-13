
# Import ELS
zip_file <- here::here("data", "ELS_2002-12_PETS_v1_0_Student_R_Datasets.zip")
load(unz(zip_file, filename = "els_02_12_byf3pststu_v1_0.rdata"))
# Note: data imported as an object "els_02_12_byf3pststu_v1_0"

els <- els_02_12_byf3pststu_v1_0 %>%
    rename(stu_id = STU_ID, sch_id = SCH_ID, sex = BYS14,
           i1 = BYS89A, i2 = BYS89B, i3 = BYS89L, i4 = BYS89R, i5 = BYS89U,
           i1_2 = F1S18A, i2_2 = F1S18B, i3_2 = F1S18C,
           i4_2 = F1S18D, i5_2 = F1S18E,
           dropout = F1EVERDO)

m_items <- c(paste0("i", 1:5), paste0("i", 1:5, "_2"))  # names of items

# helper function to recode not-appropriate answers to NA
replace_w_na_els <- function(.x, .v = c(-9, -8, -7, -6, -4, -2)) {
  .x[.x %in% .v] <- NA
  return(.x)
}
els[, m_items] <- apply(els[, m_items], MARGIN = 2, replace_w_na_els)
# Original coding: 1 = Male, 2 = Female
els$sex <- 2 - replace_w_na_els(els$sex)  # 0 = Female, 1 = Male

# Import HSLS
zip_file2 <- here::here("data", "HSLS_2017_PETS_SR_v1_0_R_Datasets.zip")
load(unz(zip_file2, filename = "hsls_17_student_pets_sr_v1_0.rdata"))
# Note: data imported as an object "hsls_17_student_pets_sr_v1_0"

hsls <- hsls_17_student_pets_sr_v1_0 %>%
    rename(stu_id = STU_ID, sch_id = SCH_ID, sex = X1SEX, 
           i1 = S1MTESTS, i2 = S1MTEXTBOOK, i4 = S1MASSEXCL, i5 = S1MSKILLS,
           i1_2 = S2MTESTS, i2_2 = S2MTEXTBOOK,
           i4_2 = S2MASSEXCL, i5_2 = S2MSKILLS,
           dropout = X2EVERDROP)

hsls$sex <- as.numeric(hsls$sex)
hsls[hsls$sex == 5, "sex"] <- NA
hsls$sex[hsls$sex == 2] <- 0 # M 1, F 0

hsls$dropout <- as.integer(hsls$dropout == "Yes")  # to be consistent with ELS

# HSLS items are coded opposite to ELS. Recode to match ELS coding.
recode_resp_hsls <- function(.x){
  .x <- as.character(.x)
  case_match(.x, "Strongly disagree" ~ "1", "Disagree" ~ "2",
             "Agree" ~ "3", "Strongly agree" ~ "4", "Missing" ~ NA, 
             "Unit non-response" ~ NA, "Item legitimate skip/NA" ~ NA, 
             .default = .x)
  }
hsls[, m_items[-c(3, 8)]] <-
  apply(hsls[, m_items[-c(3, 8)]], MARGIN = 2, recode_resp_hsls)

# create variable indicating sample
els$sample <- "ELS"
hsls$sample <- "HSLS"

# convert non-numeric variables before the merge to avoid incompatibility
char_var <- c("stu_id", "sch_id", "sex")
els[, char_var] <- apply(els[, char_var], MARGIN = 2, as.character) 
hsls[, char_var] <- apply(hsls[, char_var], MARGIN = 2, as.character) 
# convert the math items to numeric to remove attributes and labels
els[, m_items] <- apply(els[, m_items], MARGIN = 2, as.numeric) 
hsls[, m_items[-c(3, 8)]] <- apply(hsls[, m_items[-c(3, 8)]], MARGIN = 2, as.numeric)
  
# merge the ELS and HSLS data
dat_full <- full_join(els, hsls)
  
# reorder columns for ease of use
dat <- dat_full %>% select(stu_id, sch_id, sample, all_of(m_items), sex,
                           everything())
dat$sex <- as.numeric(dat$sex)

## dropout variable
# ELS: F1EVERDO # as of spring 2004
# HSLS: X2EVERDO # as of spring 2012

saveRDS(dat, 'rds/dat.rds')  # save cleaned data

# remove objects that are no longer needed from the environment
rm(dat_full, els, els_02_12_byf3pststu_v1_0, hsls, 
   hsls_17_student_pets_sr_v1_0, char_var)
