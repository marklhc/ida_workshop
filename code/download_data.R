
library(httr)
dir.create("data", showWarnings = FALSE)

# Download the ELS data (only if it doesn't already exist)
zip_file <- here("data", "ELS_2002-12_PETS_v1_0_Student_R_Datasets.zip")
if (!file.exists(zip_file)) {
    url_els <- "https://nces.ed.gov/datalab/files/zip/OnlineCodebook/ELS_2002-12_PETS_v1_0_Student_R_Datasets.zip"
    httr::GET(url_els, write_disk(zip_file, overwrite = TRUE))
}

# Download the HSLS data (only if it doesn't already exist)
zip_file2 <- here("data", "HSLS_2017_PETS_SR_v1_0_R_Datasets.zip")
if (!file.exists(zip_file2)) {
    url_hsls <- "https://nces.ed.gov/datalab/files/zip/OnlineCodebook/HSLS_2017_PETS_SR_v1_0_R_Datasets.zip"
    httr::GET(url_hsls, write_disk(zip_file2, overwrite = TRUE))
}
