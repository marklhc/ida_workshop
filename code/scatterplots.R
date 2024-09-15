
cols <- c("#F26178", "#FFCC00")

# build a ggplot scatterplot with an optional 45 degree line
# and other specifications
scatter_p <- function(dataset, x, y, xlab, ylab, title,
                      ref_line = TRUE) {
  p <- ggplot(dataset, aes(x, y, col = sample)) +
    xlab(xlab) +  ylab(ylab) + ggtitle(title) +
    theme_classic() + geom_point() + scale_color_manual(values = cols)
  if (ref_line) {
    p <- p + geom_abline(slope = 1, intercept = 0) +
      coord_fixed()
  }
  p + geom_smooth(method = "lm", linewidth = 0.5, se = FALSE)
}

#### Mean vs. factor scores

# M_B1 vs. M1
print(scatter_p(df, df$mean_score, df$partial_cont,
          ylab = "FS partial (continuous)",
          xlab = "Mean scores",
          title = "Correlation of mean scores and FS partial (cont.; M1)",
          ref_line = FALSE))

# M_B1 vs. M2
print(scatter_p(df, df$mean_score, df$approx_cont,
          ylab = "FS approximate (continuous)",
          xlab = "Mean scores",
          title = "Correlation of mean scores and FS approximate (cont.; M2)",
          ref_line = FALSE))

# M_B1 vs. M3
print(scatter_p(df, df$mean_score, df$approx_ord,
          ylab = "FS approximate (ordinal)",
          xlab = "Mean scores",
          title = "Correlation of mean scores and FS approximate (ordinal; M3)",
        ref_line = FALSE))


#### Partial vs. approximate invariance


# M1 vs. M2
print(scatter_p(df, df$partial_cont, df$partial_cont, 
          xlab = "FS partial (continuous)",
          ylab = "FS approximate (continuous)",
          "Correlation of FS: partial (cont.; M1) and approximate (cont.; M2)"))



# M4 vs. M5
print(scatter_p(score_df_noNAs, as.numeric(score_df_noNAs$partial_cont), 
          as.numeric(score_df_noNAs$approx_cont),
          xlab = "FS partial (continuous)",
          ylab = "FS approximate (continuous)",
          "[No NAs] Correlation of FS: partial (cont.; M4) and approximate (cont.; M5)"))



# M6 vs. M7
print(scatter_p(score_df_noNAs, as.numeric(score_df_noNAs$partial_ord), 
          as.numeric(score_df_noNAs$approx_ord),
          xlab = "FS partial (ordinal)",
          ylab = "FS approximate (ordinal)",
          "[No NAs] Correlation of FS: partial (ord.; M6) and approximate (ord.; M7)"))

#### Treating data as continuous vs. ordinal

# M2 vs. M3
print(scatter_p(df, df$approx_cont,
          df$approx_ord,
          xlab = "FS approximate (continuous)",
          ylab = "FS approximate (ordinal)",
          "Correlation of FS: approximate (cont.; M2) and \n approximate (ord.; M3)"))

# M4 vs. M6
print(scatter_p(score_df_noNAs, score_df_noNAs$partial_cont,
                score_df_noNAs$partial_ord,
          xlab = "FS partial (continuous)",
          ylab = "FS partial (ordinal)",
          "[No NAs] Correlation of FS: partial (cont.; M4) and \n partial (ord.; M6)"))

# M5 vs. M7
print(scatter_p(score_df_noNAs, as.numeric(score_df_noNAs$approx_cont),
          as.numeric(score_df_noNAs$approx_ord),
          xlab = "FS approximate (continuous)",
          ylab = "FS approximate (ordinal)",
          "[No NAs] Correlation of FS: approx. (cont.; M5) and \n approx. (ord.; M7)"))
