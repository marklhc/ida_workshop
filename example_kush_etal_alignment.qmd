# Load packages
library(mirt)
library(ggplot2)
library(umx)
library(R2spa)

# Data import (from https://github.com/jmk7cj/SEM-mnlfa)
data <- read.csv(here::here("SEM-mnlfa", "data.csv"))
# Sort data by Study ID
data <- data[order(data$study_id), ]
head(data)
# Sample sizes
# Study 1 = FAST; Study 2 = LIFT; Study 3 = PIRC1;
# Study 4 = PIRC2; Study 5 = SAFE
table(data[c("study_id", "race")])
# Note that the sample sizes are too small for race = black for the LIFT study,
# and for race = white for the SAFE study, so we do not include those in
# the configural model
item_names <- names(data)[5:13]

# Configural model (using `mirt`)
# 1. Create grouping variable combining study_id, sex, and race
# data$group <- paste(data$study_id, data$sex, data$race, sep = "_")
data$group <- data$study_id
group_names <- unique(data$group)
# 1b. Select only groups with enough data
selected_groups <- which(table(data$group) >= 50)
# data_sub <- subset(data, subset = group %in% names(selected_groups))
# 2. Run separate mirts (2PL)
# mirt::multipleGroup(data_sub[5:13], group = data_sub$group)
# Initialize output
# a. number of observations per group and item
item_n <- matrix(NA,
    nrow = length(selected_groups),
    ncol = length(item_names),
    dimnames = list(names(selected_groups), item_names))
# b. matrix of loadings and intercepts
lambda <- nu <- item_n
# for (g in names(selected_groups)) {
#     datg <- data[data$group == g, item_names]
#     # Drop completely missing items
#     item_n[g, ] <- colSums(!is.na(datg))
#     itemsg <- which(colSums(!is.na(datg)) > 0)
#     modg <- mirt(datg[names(itemsg)], itemtype = "2PL")
#     coefg <- coef(modg, simplify = TRUE)
#     lambda[g, names(itemsg)] <- coefg$items[, "a1"]
#     nu[g, names(itemsg)] <- coefg$items[, "d"]
# }
# The coefficients are not stable.
# We use a Bayesian prior to obtain stablized coefficients
# bmod <- "
#   f = 1-[I]
#   PRIOR = (1-[I], a1, lnorm, 0, 1), (1-[I], d, norm, 0, 5)
# "
for (g in names(selected_groups)) {
    datg <- data[data$group == g, item_names]
    # Drop completely missing items
    item_n[g, ] <- colSums(!is.na(datg))
    itemsg <- which(colSums(!is.na(datg)) > 0)
    # modg <- gsub("\\[I\\]", replacement = length(itemsg), bmod)
    fitg <- mirt(datg[names(itemsg)], model = 1,
                 itemtype = "2PL")
    coefg <- coef(fitg, simplify = TRUE)
    lambda[g, names(itemsg)] <- coefg$items[, "a1"]
    nu[g, names(itemsg)] <- coefg$items[, "d"]
}
# 3. Alignment optimization (not converging)
aligned <- sirt::invariance.alignment(
    lambda = lambda,
    nu = nu,
    wgt = sqrt(item_n),
    optimizer = "nlminb",
    eps = .0001,
    control = list(maxit = 10000)
)
# 4. Effect size
thres <- aligned$nu.aligned[, c(1, 3:9)]
colnames(thres) <- seq_len(ncol(thres))
# pooled SD
vars <- apply(data[item_names],
    MARGIN = 2,
    FUN = \(x) tapply(x,
        INDEX = data$study_id,
        FUN = var, na.rm = TRUE
    )
)
item_n2 <- ifelse(item_n > 0, item_n, NA)
pooled_sd <- sqrt(
    colSums(vars * (item_n2 - 1), na.rm = TRUE) /
        colSums(item_n2 - 1, na.rm = TRUE)
)
# Noninvariance effect size for items 1, 3-9
pinsearch::fmacs_ordered(
    thresholds = thres,
    loadings = aligned$lambda.aligned[, c(1, 3:9)],
    link = "logit",
    pooled_item_sd = pooled_sd[c(1, 3:9)]
)
# Noninvariance effect size for item 2
thres2 <- aligned$nu.aligned[-1, 2, drop = FALSE]
colnames(thres2) <- 1
pinsearch::fmacs_ordered(
    thresholds = thres2,
    loadings = aligned$lambda.aligned[-1, 2, drop = FALSE],
    link = "logit",
    pooled_item_sd = pooled_sd[2]
)
# The magnitude of noninvariance effect size is not trivial (with fmacs > .1)
# for all items.

# 5. Compute factor scores
# 5a. Refit the configural model, with latent means and variances from the
#     aligned parameters
# 5b. Save factor scores, using the same prior mean and variance
(aligned_pars <- aligned$pars)
mod <- "
  f = 1-[I]
  START = (GROUP, MEAN_1, [mu]), (GROUP, COV_11, [sigma2])
"
# Substitute into the aligned means and variances
fsg <- NULL
for (g in names(selected_groups)) {
    datg <- data[data$group == g, item_names]
    # Drop completely missing items
    itemsg <- which(item_n[g, ] > 0)
    modg <- gsub("\\[I\\]", replacement = length(itemsg), mod)
    modg <- gsub("\\[mu\\]", replacement = aligned_pars[g, 1], modg)
    modg <- gsub("\\[sigma2\\]", replacement = aligned_pars[g, 2], modg)
    fitg <- mirt(datg[names(itemsg)], model = modg,
                 itemtype = "2PL")
    fsg <- rbind(
        fsg,
        data.frame(
            fscores(fitg,
                method = "EAP", mean = 0, cov = 1, full.scores.SE = TRUE
            ),
            group = g
        )
    )
}
# Average reliability
var(fsg[, 1], na.rm = TRUE) /
    (var(fsg[, 1], na.rm = TRUE) + mean(fsg[, 2]^2, na.rm = TRUE))
# Group-specific reliability
relg <- tapply(fsg, INDEX = fsg$group,
    FUN = \(x) var(x[, 1], na.rm = TRUE) /
        (var(x[, 1], na.rm = TRUE) + mean(x[, 2]^2, na.rm = TRUE)))

# 6. Naive second-stage analysis
# 6a. Descriptives
ggplot(fsg, aes(x = f)) +
    geom_histogram() +
    facet_wrap(~ group) +
    labs(x = "EAP score")
# 6b. Merged with other variables
# Warning: the code requires data to be sorted by study_id before analyses
data <- cbind(data, fsg)
# 6c. Recode discrete variable to be used with OpenMx
data$hs2 <- mxFactor(data$hs, levels = c(0, 1))
# 6d. Pooled analysis without adjusting for unreliability
# TOCA-R -> high school completion
probitreg_umx <- umxRAM("
  hs2 ~ f
  hs2 ~ 0 * 1
  hs2 ~~ 1 * hs2
",
    data = data)
plot(probitreg_umx)  # path diagram

# 7. Second stage analysis, adjusting for measurement error
# 7b. Define loading and error variances
# Loading = Latent variance (for scoring) - SE^2
data$rel_f <- 1 - data$SE_f^2
data$ev_f <- data$SE_f^2 * data$rel_f
# 7c. Drop rows with missing rel_f
data2 <- data[!is.na(data$rel_f), ]
# 7d. Incorporate measurement model for factor scores
# Note: as EAP scores are only used as a predictor, the
#       unstandardized coefficients are not affected by
#       the measurement model.
probitreg_2spa <- umxRAM(
    "2spa",
    umxPath("eta", to = "hs2", values = -0.23),
    umxPath("eta", to = "f", labels = "data.rel_f", free = FALSE),
    umxPath(var = "f", labels = "data.ev_f", free = FALSE),
    umxPath(v.m. = "eta", values = c(1.317, -0.784)),
    umxPath(v1m0 = "hs2"),
    data = data2,
    tryHard = "ordinal"
)

# 8. Adding grouping variable
data$gp <- as.numeric(factor(data$group))
probitreg_umx <- umxRAM(
    "fspa",
    umxPath("f", to = "hs2"),
    umxPath(v1m0 = "hs2"),
    umxPath(v.m. = "f"),
    data = data,
    tryHard = "ordinal",
    group = "gp")

# Optimization by hand
# sim_fun <- function(x, lam, nu, wgt, eps = 0.001) {
#     n_gps <- nrow(lam)
#     dpsi <- x[seq_len(n_gps - 1)]
#     dalpha <- x[seq_len(n_gps - 1) + n_gps - 1]
#     new_lam <- lam / (1 + c(0, dpsi))
#     new_nu <- nu - c(0, dalpha) * new_lam
#     tnew_lam <- t(new_lam)
#     tnew_nu <- t(new_nu)
#     loss <- 0
#     for (g in seq_len(n_gps - 1)) {
#         g1 <- (g + 1):n_gps
#         loss <- loss +
#             sum(
#                 sqrt(sqrt((tnew_lam[, g1] - tnew_lam[, g])^2 + eps)),
#                 sqrt(sqrt((tnew_nu[, g1] - tnew_nu[, g])^2 + eps)),
#                 na.rm = TRUE
#             )
#     }
#     loss
# }

# opt <- optim(rep(0.01, 30), sim_fun,
#              lam = lambda, nu = nu, wgt = sqrt(item_n), eps = 0.0001,
#              method = "L-BFGS-B",
#              control = list(maxit = 10000))

# Exclude Study 1
sirt::invariance.alignment(
    lambda = lambda[-(1:4), ],
    nu = nu[-(1:4), ],
    wgt = sqrt(item_n[-(1:4), ]),
    optimizer = "nlminb",
    control = list(iter.max = 1000)
)

library(blavaan)
bcfa("
  f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
", data = data[data$group %in% names(selected_groups), item_names] |>
    setNames(c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")),
    ordered = TRUE,
    n.chains = 1, burnin = 200, sample = 200)