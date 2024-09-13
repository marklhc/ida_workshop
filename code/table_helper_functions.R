# helper function to build a table of mean and variances of scores by approach and dataset
mnsd <- function(vals) {
  c(mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE))
}

build_mn_sd_tab <- function(scores, scores_noNAs) {
mn_sd_tab <- as.data.frame(rbind(c("M1", "", "M2", "", "", "", "M3", ""),
  cbind(
    round(rbind(c(mnsd(scores[scores$sample == "ELS", "partial_cont"]),
                  mnsd(scores[scores$sample == "ELS", "approx_cont"])),
                c(mnsd(scores[scores$sample == "HSLS", "partial_cont"]),
                  mnsd(scores[scores$sample == "HSLS", "approx_cont"])),
                c(mnsd(scores$partial_cont), mnsd(scores$approx_cont))), 3),
    cbind(rep("-", 3), rep("-", 3)), 
    round(rbind(mnsd(scores[scores$sample == "ELS", "approx_ord"]),
                mnsd(scores[scores$sample == "HSLS", "approx_ord"]),
                mnsd(scores$approx_ord)), 3)),
  c("M4", "", "M5", "", "M6", "", "M7", ""),
  round(rbind(
    c(mnsd(scores_noNAs[scores_noNAs$sample == "ELS", "partial_cont"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "ELS", "approx_cont"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "ELS", "partial_ord"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "ELS", "approx_ord"])),
    c(mnsd(scores_noNAs[scores_noNAs$sample == "HSLS", "partial_cont"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "HSLS", "approx_cont"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "HSLS", "partial_ord"]),
      mnsd(scores_noNAs[scores_noNAs$sample == "HSLS", "approx_ord"])),
    c(mnsd(scores_noNAs$partial_cont), 
      mnsd(scores_noNAs$approx_cont), 
      mnsd(scores_noNAs$partial_ord), 
      mnsd(scores_noNAs$approx_ord))), 3)
  ))

mean_sc_mnsd <- 
  rbind(c("M_B1", ""),
        cbind(round(rbind(mnsd(scores[scores$sample=="ELS", "mean_score"]), 
                          mnsd(scores[scores$sample=="HSLS", "mean_score"]),
                          mnsd(scores$mean_score)), 3)),
        c("M_B1", ""),
        cbind(round(rbind(mnsd(scores_noNAs[scores_noNAs$sample=="ELS", "mean_score"]), 
                          mnsd(scores_noNAs[scores_noNAs$sample=="HSLS", "mean_score"]),
                          mnsd(scores_noNAs$mean_score)), 3)))
mn_sd_tab <- cbind(mean_sc_mnsd, mn_sd_tab)

rownames(mn_sd_tab) <- c("N = 30,740", "ELS", "HSLS", "overall", "N = 29,202 ", 
                         "ELS ", "HSLS ", "overall ")
colnames(mn_sd_tab) <- rep(c("M", "SD"), 5)

mn_sd_tab_kbl <-
  kbl(mn_sd_tab, align = "c", digits = 3, booktabs = TRUE,
    caption = "Mean and SD of scores by approach and dataset") %>% kable_classic() %>%
  add_header_above(c("", "Mean Score" = 2, "Partial" = 2, "Approximate" = 2,
                     "Partial" = 2, "Approximate" = 2)) %>%
  add_header_above(c(" " = 3, "Continuous" = 4, "Ordinal" = 4)) %>%
  row_spec(5, bold = T, extra_css = "border-bottom: 1px solid") %>% 
  row_spec(4, extra_css = "border-bottom: 1px solid") %>% 
  row_spec(1, bold = T, extra_css = "border-bottom: 1px solid") 
return(mn_sd_tab_kbl)
}


build_loadings_tab <- function() {
  loadings_tab <- 
    cbind(est_partial$ELS$lambda, 
          c(est_partial$HSLS$lambda[1:2], NA, est_partial$HSLS$lambda[3:4]),
          est_align$ELS$lambda, 
          c(est_align$HSLS$lambda[1:2], NA, est_align$HSLS$lambda[3:4]),
          est_partial_ord$ELS$lambda, c(est_partial_ord$HSLS$lambda[1:2], NA,
                                        est_partial_ord$HSLS$lambda[3:4]),
          est_align_ord[est_align_ord$group=="ELS", "a1"],
          est_align_ord[est_align_ord$group=="HSLS", "a1"])
  colnames(loadings_tab) <- NULL 
  
  loadings_tab_kbl <- 
    kbl(loadings_tab, booktabs = T, align = "c", table.attr = "style='width:60%;'",
      caption = "Loading estimates by approach and dataset", longtable = TRUE, digits = 2) %>%
    kable_classic() %>% 
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1,
                       "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c("", "Partial (M1)" = 2, "Approximate (M2)" = 2,
                       "Partial (M6)" = 2, "Approximate (M3)" = 2)) %>%
    add_header_above(c("  " = 1, "Continuous" = 4, "Ordinal" = 4))
  return(loadings_tab_kbl)
}
  
build_intercepts_tab <- function() {
  intercepts_tab <- 
    cbind(est_partial$ELS$nu, 
          c(est_partial$HSLS$nu[1:2], NA, est_partial$HSLS$nu[3:4]),
          est_align$ELS$nu, 
          c(est_align$HSLS$nu[1:2], NA, est_align$HSLS$nu[3:4]))
  
  colnames(intercepts_tab) <- NULL 
  
  intercepts_tab_kbl <- kbl(intercepts_tab, booktabs = T, align = "c", table.attr = "style='width:60%;'",
      caption = "Intercept estimates by approach and dataset", longtable = TRUE, digits = 2) %>%
    kable_classic() %>% 
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c("", "Partial (M1)" = 2, "Approximate (M2)" = 2)) %>%
    add_header_above(c("  " = 1, "Continuous" = 4))
  return(intercepts_tab_kbl)  
}


build_thresholds_tab <- function() {
  thresholds_tab <- cbind(est_partial_ord$ELS$tau, 
                         c(est_partial_ord$HSLS$tau[1:6], rep(NA, 3),
                           est_partial_ord$HSLS$tau[7:12]), t(est_align_ord_th))
  colnames(thresholds_tab) <- NULL 
  
  thresholds_tab_kbl <- kbl(thresholds_tab, booktabs = T, align = "c", table.attr = "style='width:60%;'",
      caption = "Threshold estimates by approach and dataset", longtable = TRUE, digits = 2) %>%
    kable_classic() %>% 
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c("", "Partial (M6)" = 2, "Approximate (M3)" = 2)) %>%
    add_header_above(c("  " = 1, "Ordinal" = 4)) %>%
    row_spec(c(3, 6, 9, 12), extra_css = "border-bottom: 1px solid")
  
  return(thresholds_tab_kbl)
}


build_latent_tab <- function() {
  alpha_psi_tab <- rbind(
    cbind(est_partial$ELS$alpha, est_partial$HSLS$alpha,
          est_align$ELS$alpha, est_align$HSLS$alpha,
          est_partial_ord$ELS$alpha, est_partial_ord$HSLS$alpha,
          mirt.wrapper.coef(mod_aligned)$GroupPars$ELS[1], 
          mirt.wrapper.coef(mod_aligned)$GroupPars$HSLS[1]
    ),
    cbind(est_partial$ELS$psi, est_partial$HSLS$psi,
          est_align$ELS$psi, est_align$HSLS$psi,
          est_partial_ord$ELS$psi, est_partial_ord$HSLS$psi,
          (mirt.wrapper.coef(mod_aligned)$GroupPars$ELS[2])^2, 
          (mirt.wrapper.coef(mod_aligned)$GroupPars$HSLS[2])^2
    ))
  colnames(alpha_psi_tab) <- NULL
  rownames(alpha_psi_tab) <- c("Latent mean", "Latent variance")
  latent_tab_kbl <- 
    kbl(alpha_psi_tab, booktabs = T, align = "c", table.attr = "style='width:60%;'",
      caption = "Latent mean and variance estimates by approach and dataset", longtable = TRUE, digits = 2) %>%
    kable_classic() %>% 
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1,
                       "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c("", "Partial (M1)" = 2, "Approximate (M2)" = 2,
                       "Partial (M6)" = 2, "Approximate (M3)" = 2)) %>%
    add_header_above(c("  " = 1, "Continuous" = 4, "Ordinal" = 4))
  return(latent_tab_kbl)  
}
