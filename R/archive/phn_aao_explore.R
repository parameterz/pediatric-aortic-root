# =====================================================================
#  PHN ascending aorta (AAo) exploration
#   1. How does AAo diameter scale with height (and BSA)?
#   2. Candidate z-score models built on the diameter itself
#   3. Are those z-scores independent of body size and age?
#   4. Stress test of CSA/height (Chaudhary et al., Pediatr Cardiol 2026)
#   5. Agreement, error sensitivity, and z = 2 cutoffs across sizes
#
#  Expects the data frame `echo_zscores_phn` with:
#    AR_AAD (cm), ECHO_AGE (y), ECHO_HT (cm), ECHO_WT (kg), BSA_Haycock (m2)
#  plus a sex column for the Chaudhary sections (see Section 0).
#
#  Packages: ggplot2, mgcv (mgcv ships with R)
#  Output: console + aao_explore/report.txt + PNG figures
# =====================================================================

library(ggplot2)
library(mgcv)

out_dir <- "aao_explore"
dir.create(out_dir, showWarnings = FALSE)
sink(file.path(out_dir, "report.txt"), split = TRUE)   # sink() at the end closes it

show <- function(title, x = NULL) {
  cat("\n==== ", title, " ====\n", sep = "")
  if (!is.null(x)) print(x)
}
save_plot <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  print(p)
}
skew <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }


# ---------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------
raw <- echo_zscores_phn

# Sex is needed for the Chaudhary model. If auto-detection fails, set
# sex_col by hand. CHECK THE CODING in the printed table: numeric
# coding is ASSUMED to be 1 = male.
sex_col <- intersect(c("SEX", "Sex", "sex", "GENDER", "Gender", "gender",
                       "PAT_SEX", "ECHO_SEX"), names(raw))[1]

d <- data.frame(
  aad = raw$AR_AAD,            # cm
  age = raw$ECHO_AGE,          # years
  ht  = raw$ECHO_HT / 100,     # m
  wt  = raw$ECHO_WT,           # kg
  bsa = raw$BSA_Haycock        # m2
)

has_sex <- !is.na(sex_col)
if (has_sex) {
  s <- raw[[sex_col]]
  show(paste("Sex coding in column", sex_col), table(s, useNA = "ifany"))
  d$male <- if (is.numeric(s)) as.integer(s == 1) else
                               as.integer(grepl("^m", tolower(as.character(s))))
} else {
  message("No sex column found -- set sex_col by hand. Chaudhary sections need it.")
  d$male <- NA_integer_
}

d <- d[complete.cases(d[, c("aad", "age", "ht", "wt", "bsa")]), ]
if (has_sex) d <- d[!is.na(d$male), ]
d$sexf <- factor(d$male, 0:1, c("F", "M"))
d$bmi  <- d$wt / d$ht^2
d$lad  <- log(d$aad)
d$lht  <- log(d$ht)
d$lbsa <- log(d$bsa)
d$ht_decile <- cut(d$ht, quantile(d$ht, 0:10 / 10), include.lowest = TRUE, labels = 1:10)
# Chaudhary's age bins (min <= age < max); 18.0 goes in the last bin
d$age_grp <- cut(d$age, c(-Inf, 1/12, 3, 6, 12, 16, Inf), right = FALSE,
                 labels = c("<1mo", "1mo-3y", "3-6y", "6-12y", "12-16y", "16-18y"))

show(paste("n =", nrow(d)), summary(d[, c("aad", "age", "ht", "wt", "bsa", "bmi")]))
show("Correlations among size variables (collinearity check)",
     round(cor(d[, c("lht", "age", "lbsa", "wt")]), 3))


# ---------------------------------------------------------------------
# 1. Scaling of AAo diameter with height and BSA
# ---------------------------------------------------------------------
fit_exp <- function(f, data) {
  m <- lm(f, data)
  ci <- confint(m)[2, ]
  c(n = nrow(data), exponent = unname(coef(m)[2]), lo95 = ci[[1]], hi95 = ci[[2]],
    r2 = summary(m)$r.squared, resid_sd_log = sigma(m))
}
exps <- rbind(
  `height, all`     = fit_exp(lad ~ lht,  d),
  `height, <2 y`    = fit_exp(lad ~ lht,  subset(d, age < 2)),
  `height, >=2 y`   = fit_exp(lad ~ lht,  subset(d, age >= 2)),
  `height, >=12 y`  = fit_exp(lad ~ lht,  subset(d, age >= 12)),
  `BSA, all`        = fit_exp(lad ~ lbsa, d),
  `BSA, >=2 y`      = fit_exp(lad ~ lbsa, subset(d, age >= 2))
)
show("1a. Allometric exponents (log AAD ~ log size)", round(exps, 3))

# Growth vs static allometry: the height exponent WITHIN a narrow age band
# (kids of the same age, differing in size) vs ACROSS ages (growth).
d$age_yr <- factor(pmin(floor(d$age), 17))
within_age <- coef(lm(lad ~ lht + age_yr, d))["lht"]
by_band <- t(sapply(split(d, d$age_grp), function(s) fit_exp(lad ~ lht, s)))
show("1b. Height exponent within 1-year age bands (pooled)", round(within_age, 3))
show("1c. Height exponent within each Chaudhary age group", round(by_band, 3))

# Local exponent along the height range (derivative of a smooth on log-log)
g_ht <- gam(lad ~ s(lht, k = 10), data = d, method = "REML")
grid <- data.frame(lht = seq(min(d$lht), max(d$lht), length.out = 200))
eps <- 1e-4
Xd <- (predict(g_ht, transform(grid, lht = lht + eps), type = "lpmatrix") -
       predict(g_ht, transform(grid, lht = lht - eps), type = "lpmatrix")) / (2 * eps)
grid$slope <- drop(Xd %*% coef(g_ht))
grid$se    <- sqrt(rowSums((Xd %*% vcov(g_ht)) * Xd))
grid$ht_cm <- exp(grid$lht) * 100

save_plot(
  ggplot(grid, aes(ht_cm, slope)) +
    geom_ribbon(aes(ymin = slope - 2 * se, ymax = slope + 2 * se), alpha = 0.2) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1,   linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    labs(x = "Height (cm)", y = "Local exponent: d log(AAD) / d log(Ht)",
         title = "Local height exponent of AAo diameter",
         subtitle = "dashed = 1 (geometric, growth); dotted = 0.5 (adult-like)"),
  "01_local_exponent")

save_plot(
  ggplot(d, aes(ht * 100, aad * 10)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    geom_smooth(method = "gam", formula = y ~ s(x), colour = "firebrick", se = FALSE) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Height (cm, log)", y = "AAo diameter (mm, log)",
         title = "AAo diameter vs height, log-log",
         subtitle = "black = power law; red = smooth"),
  "02_loglog_scatter")

# Does age (or sex) add anything once height is in?
m0 <- lm(lad ~ lht, d)
m1 <- lm(lad ~ lht + age, d)
show("1d. Does age add to log height?", anova(m0, m1))
show("    coefficients", round(summary(m1)$coefficients, 4))
cat(sprintf("    age effect at fixed height: %+.2f%% per 10 years\n",
            100 * (exp(10 * coef(m1)["age"]) - 1)))
if (has_sex) {
  m2 <- lm(lad ~ lht + age + male, d)
  show("1e. Does sex add to log height + age?", anova(m1, m2))
  cat(sprintf("    male effect at fixed height/age: %+.2f%%\n",
              100 * (exp(coef(m2)["male"]) - 1)))
}


# ---------------------------------------------------------------------
# 2. Candidate z-score models on the diameter
# ---------------------------------------------------------------------
# M1: power law, constant SD on the log scale  (== AAD / Ht^b)
m_pow <- lm(lad ~ lht, d)
a_pow <- coef(m_pow)[[1]]; b_pow <- coef(m_pow)[[2]]; s_pow <- sigma(m_pow)
d$z_pow <- resid(m_pow) / s_pow
show("2a. Power-law model (log scale)",
     round(c(intercept = a_pow, exponent = b_pow, sd_log = s_pow), 4))

# Heteroscedasticity check for M1
het <- lm(abs(resid(m_pow)) ~ d$lht)
show("2b. |residual| vs log height (slope ~ 0 means constant log-SD)",
     round(summary(het)$coefficients, 4))

# M2: smooth mean AND smooth SD on log height (location-scale GAM)
m_ls <- gam(list(lad ~ s(lht, k = 10), ~ s(lht, k = 6)),
            family = gaulss(), data = d, method = "REML")
ls_pred <- function(newdata) {
  p <- predict(m_ls, newdata = newdata, type = "response")
  sd <- if (median(p[, 2]) > 1) 1 / p[, 2] else p[, 2]   # col 2 may be 1/sd
  data.frame(mu = p[, 1], sd = sd)
}
lp <- ls_pred(d)
d$z_ls <- (d$lad - lp$mu) / lp$sd
show("2c. Location-scale GAM: range of fitted log-SD", round(range(lp$sd), 4))

# M3: Mahgerefteh 2021 published height model (AAD/Ht^1, kids >= 2 y)
d$z_mah <- (d$aad / d$ht - 1.39) / 0.14


# ---------------------------------------------------------------------
# 3. Chaudhary et al. CSA/ht model (values copied from their supplement)
# ---------------------------------------------------------------------
d$csa_ht  <- pi * (d$aad / 2)^2 / d$ht          # cm2/m
d$csa_ht2 <- pi * (d$aad / 2)^2 / d$ht^2        # cm2/m2

if (has_sex) {
  ch_bins <- data.frame(
    age_grp = levels(d$age_grp),
    mean_F = c(1.1566680, 1.3894141, 1.6674859, 1.9498509, 2.3370879, 2.4963250),
    sd_F   = c(0.2718176, 0.3164411, 0.3513082, 0.3534838, 0.4802758, 0.5387956),
    mean_M = c(1.2500483, 1.3956703, 1.7176683, 2.0984815, 2.3480720, 2.6417850),
    sd_M   = c(0.2862416, 0.3201912, 0.3666948, 0.4530524, 0.4751761, 0.5027709))
  ch_sheet <- c(int = 1.2612332551, age = 0.0538404716,
                male = 0.0751494789, male_age = 0.0031328790)
  ch_paper <- replace(ch_sheet, "age", 0.08)      # abstract reports beta(age) = 0.08
  ch_mu <- function(cf, age, male)
    unname(cf["int"] + cf["age"] * age + cf["male"] * male + cf["male_age"] * age * male)

  i <- match(as.character(d$age_grp), ch_bins$age_grp)
  d$ch_bin_mean <- ifelse(d$male == 1, ch_bins$mean_M[i], ch_bins$mean_F[i])
  d$ch_bin_sd   <- ifelse(d$male == 1, ch_bins$sd_M[i],   ch_bins$sd_F[i])
  d$z_ch_cat    <- (d$csa_ht - d$ch_bin_mean) / d$ch_bin_sd
  d$z_ch_sheet  <- (d$csa_ht - ch_mu(ch_sheet, d$age, d$male)) / d$ch_bin_sd
  d$z_ch_paper  <- (d$csa_ht - ch_mu(ch_paper, d$age, d$male)) / d$ch_bin_sd

  # 3a. Can we reproduce their bin table from this data?
  rep_tab <- do.call(rbind, lapply(split(d, list(d$sexf, d$age_grp)), function(s)
    data.frame(sex = s$sexf[1], age_grp = s$age_grp[1], n = nrow(s),
               mean = mean(s$csa_ht), sd = sd(s$csa_ht))))
  j <- match(as.character(rep_tab$age_grp), ch_bins$age_grp)
  rep_tab$their_mean <- ifelse(rep_tab$sex == "M", ch_bins$mean_M[j], ch_bins$mean_F[j])
  rep_tab$their_sd   <- ifelse(rep_tab$sex == "M", ch_bins$sd_M[j],   ch_bins$sd_F[j])
  show("3a. Their bin means/SDs vs recomputed from this data",
       transform(rep_tab, mean = round(mean, 3), sd = round(sd, 3),
                 their_mean = round(their_mean, 3), their_sd = round(their_sd, 3)))

  # 3b. Refit their regression: which coefficients come out?
  m_ch <- lm(csa_ht ~ age * male, d)
  show("3b. Refit csa_ht ~ age * male  (compare: sheet / paper)",
       rbind(refit = round(coef(m_ch), 4),
             sheet = round(ch_sheet, 4), paper = round(ch_paper, 4)))

  # 3c. Is age a surrogate for the missing height exponent?
  m_a  <- lm(csa_ht ~ age * male, d)
  m_h  <- lm(csa_ht ~ ht, d)
  m_ah <- lm(csa_ht ~ ht + age * male, d)
  show("3c. R2: age*sex vs height vs both",
       round(c(age_sex = summary(m_a)$r.squared, height = summary(m_h)$r.squared,
               both = summary(m_ah)$r.squared), 3))
  show("    age/sex terms once height is in the model",
       round(summary(m_ah)$coefficients, 4))
  show("    exponent of CSA/ht on height (0 = size-neutral)",
       round(fit_exp(log(csa_ht) ~ lht, d), 3))
  show("    exponent of CSA/ht^2 on height (0 = size-neutral)",
       round(fit_exp(log(csa_ht2) ~ lht, d), 3))
  show("    csa/ht^2 ~ age * male (do age/sex effects vanish?)",
       round(summary(lm(csa_ht2 ~ age * male, d))$coefficients, 4))

  # Their mean lines vs the data
  age_grid <- expand.grid(age = seq(0, 18, by = 0.1), male = 0:1)
  age_grid$sexf <- factor(age_grid$male, 0:1, c("F", "M"))
  age_grid$sheet <- ch_mu(ch_sheet, age_grid$age, age_grid$male)
  age_grid$paper <- ch_mu(ch_paper, age_grid$age, age_grid$male)
  save_plot(
    ggplot(d, aes(age, csa_ht)) +
      geom_point(alpha = 0.2, size = 0.7) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 10), colour = "black", se = FALSE) +
      geom_line(data = age_grid, aes(y = sheet), colour = "firebrick", linewidth = 1) +
      geom_line(data = age_grid, aes(y = paper), colour = "steelblue", linewidth = 1,
                linetype = "dashed") +
      facet_wrap(~ sexf) +
      labs(x = "Age (y)", y = "CSA/ht (cm2/m)", title = "CSA/ht vs age",
           subtitle = "black = data smooth; red = spreadsheet line; blue dashed = paper beta"),
    "03_csa_ht_vs_age", w = 9)

  save_plot(
    ggplot(d, aes(ht * 100, csa_ht, colour = sexf)) +
      geom_point(alpha = 0.25, size = 0.7) +
      geom_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
      labs(x = "Height (cm)", y = "CSA/ht (cm2/m)", colour = "Sex",
           title = "CSA/ht still rises with height (not size-neutral)"),
    "04_csa_ht_vs_height")
}


# ---------------------------------------------------------------------
# 4. Are the z-scores independent of body size and age?
# ---------------------------------------------------------------------
check_z <- function(z, name, data = d) {
  covs <- c(lht = "log height", age = "age", wt = "weight", bsa = "BSA", bmi = "BMI")
  if (has_sex) covs <- c(covs, male = "male sex")
  slopes <- t(sapply(names(covs), function(v) {
    x <- data[[v]]
    xs <- if (v == "male") x else (x - mean(x)) / sd(x)   # per 1 SD of covariate
    m <- lm(z ~ xs)
    c(slope_per_SD = unname(coef(m)[2]), p = summary(m)$coefficients[2, 4],
      r2 = summary(m)$r.squared)
  }))
  rownames(slopes) <- covs
  by_tab <- function(g) t(sapply(split(z, g), function(v)
    c(n = length(v), mean = mean(v), sd = sd(v),
      pct_gt2 = 100 * mean(v > 2), pct_lt_m2 = 100 * mean(v < -2))))
  show(paste0("4. ", name, ": overall"),
       round(c(mean = mean(z), sd = sd(z), skew = skew(z),
               pct_gt2 = 100 * mean(z > 2), pct_lt_m2 = 100 * mean(z < -2)), 3))
  show(paste0("   ", name, ": z ~ covariate (slope per SD; ideal 0)"), signif(slopes, 3))
  show(paste0("   ", name, ": by height decile (ideal mean 0, sd 1, ~2.3% each tail)"),
       round(by_tab(data$ht_decile), 2))
  show(paste0("   ", name, ": by age group"), round(by_tab(data$age_grp), 2))
  invisible(slopes)
}

z_models <- c(z_pow = "Power law (log AAD ~ log Ht)",
              z_ls  = "Location-scale GAM on log Ht")
if (has_sex) z_models <- c(z_models,
              z_ch_cat   = "CSA/ht, categorical age",
              z_ch_sheet = "CSA/ht, continuous age (spreadsheet)",
              z_ch_paper = "CSA/ht, continuous age (paper beta)")
for (v in names(z_models)) check_z(d[[v]], z_models[[v]])
# Mahgerefteh is published for >= 2 y only
check_z(d$z_mah[d$age >= 2], "Mahgerefteh AAD/Ht (>= 2 y only)", d[d$age >= 2, ])

zl <- do.call(rbind, lapply(names(z_models), function(v) {
  rbind(data.frame(model = z_models[[v]], var = "Height (cm)", x = d$ht * 100, z = d[[v]]),
        data.frame(model = z_models[[v]], var = "Age (y)",     x = d$age,      z = d[[v]]),
        data.frame(model = z_models[[v]], var = "Weight (kg)", x = d$wt,       z = d[[v]]),
        data.frame(model = z_models[[v]], var = "BMI",         x = d$bmi,      z = d[[v]]))
}))
save_plot(
  ggplot(zl, aes(x, z)) +
    geom_point(alpha = 0.08, size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = "grey40") +
    geom_smooth(method = "gam", formula = y ~ s(x, k = 8), colour = "firebrick", se = FALSE) +
    facet_grid(model ~ var, scales = "free_x") +
    coord_cartesian(ylim = c(-4, 4)) +
    labs(x = NULL, y = "z", title = "z-score residual checks (red should be flat at 0)"),
  "05_z_residual_checks", w = 12, h = 2.4 * length(z_models) + 1)


# ---------------------------------------------------------------------
# 5. Agreement, error sensitivity, the 8-year-old, and z = 2 cutoffs
# ---------------------------------------------------------------------
if (has_sex) {
  show("5a. Classification agreement: power law z > 2 vs CSA/ht categorical z > 2",
       table(power_law = d$z_pow > 2, csa_ht_cat = d$z_ch_cat > 2))
  show("    power law z > 2 vs CSA/ht spreadsheet z > 2",
       table(power_law = d$z_pow > 2, csa_ht_sheet = d$z_ch_sheet > 2))

  save_plot(
    ggplot(d, aes(ht * 100, z_ch_cat - z_pow)) +
      geom_point(alpha = 0.2, size = 0.7) +
      geom_smooth(method = "gam", formula = y ~ s(x), colour = "firebrick") +
      geom_hline(yintercept = 0) +
      labs(x = "Height (cm)", y = "z (CSA/ht categorical) - z (power law)",
           title = "Where do the two z-scores disagree?"),
    "06_z_difference_vs_height")

  # Sensitivity: how much does z move for a +1 mm diameter error?
  d$dz_pow <- (log(d$aad + 0.1) - d$lad) / s_pow
  d$dz_ch  <- (pi * ((d$aad + 0.1) / 2)^2 - pi * (d$aad / 2)^2) / d$ht / d$ch_bin_sd
  sens <- t(sapply(split(d, d$age_grp), function(s)
    c(n = nrow(s), median_AAD_mm = median(s$aad) * 10,
      dz_power_law = mean(s$dz_pow), dz_csa_ht = mean(s$dz_ch))))
  show("5b. Change in z for a +1 mm diameter error, by age group", round(sens, 2))
}

# Dan's example: 8-year-old girl, 128 cm, AAo 24 mm
# (Lauzon-Schnitka z = +2.81; PHN BSA z = +3.41)
case <- data.frame(aad = 2.4, ht = 1.28, age = 8, male = 0)
case_ls <- ls_pred(data.frame(lht = log(case$ht)))
case_z <- c(
  power_law   = (log(case$aad) - a_pow - b_pow * log(case$ht)) / s_pow,
  loc_scale   = (log(case$aad) - case_ls$mu) / case_ls$sd,
  mahgerefteh = (case$aad / case$ht - 1.39) / 0.14)
if (has_sex) {
  csa <- pi * (case$aad / 2)^2 / case$ht
  case_z <- c(case_z,
    csa_ht_cat   = (csa - ch_bins$mean_F[4]) / ch_bins$sd_F[4],
    csa_ht_sheet = (csa - ch_mu(ch_sheet, 8, 0)) / ch_bins$sd_F[4],
    csa_ht_paper = (csa - ch_mu(ch_paper, 8, 0)) / ch_bins$sd_F[4])
}
show("5c. 8 y girl, 128 cm, AAo 24 mm (L-S +2.81, PHN BSA +3.41)", round(case_z, 2))

# z = 2 cutoff diameter (mm) across typical kids
if (has_sex) {
  ht_fit <- gam(ht ~ sexf + s(age, by = sexf, k = 10), data = d)
  cut_grid <- expand.grid(age = c(0.03, 0.5, 1, 2, 4, 8, 12, 15, 17), sexf = c("F", "M"))
  cut_grid$sexf <- factor(cut_grid$sexf, levels = c("F", "M"))
  cut_grid$male <- as.integer(cut_grid$sexf == "M")
  cut_grid$ht_cm <- predict(ht_fit, cut_grid) * 100
  h <- cut_grid$ht_cm / 100
  gp <- ls_pred(data.frame(lht = log(h)))
  gi <- match(as.character(cut(cut_grid$age, c(-Inf, 1/12, 3, 6, 12, 16, Inf), right = FALSE,
              labels = levels(d$age_grp))), ch_bins$age_grp)
  g_sd   <- ifelse(cut_grid$male == 1, ch_bins$sd_M[gi],   ch_bins$sd_F[gi])
  g_mean <- ifelse(cut_grid$male == 1, ch_bins$mean_M[gi], ch_bins$mean_F[gi])
  csa_to_mm <- function(csa_ht, h) 2 * sqrt(csa_ht * h / pi) * 10
  cut_grid$power_law <- exp(a_pow + b_pow * log(h) + 2 * s_pow) * 10
  cut_grid$loc_scale <- exp(gp$mu + 2 * gp$sd) * 10
  cut_grid$mahgerefteh <- ifelse(cut_grid$age >= 2, (1.39 + 2 * 0.14) * h * 10, NA)
  cut_grid$csa_cat   <- csa_to_mm(g_mean + 2 * g_sd, h)
  cut_grid$csa_sheet <- csa_to_mm(ch_mu(ch_sheet, cut_grid$age, cut_grid$male) + 2 * g_sd, h)
  cut_grid$csa_paper <- csa_to_mm(ch_mu(ch_paper, cut_grid$age, cut_grid$male) + 2 * g_sd, h)
  num <- sapply(cut_grid, is.numeric)
  cut_grid[num] <- lapply(cut_grid[num], round, 1)
  show("5d. Diameter (mm) at z = +2 for a median-height child, by method", cut_grid)
}

sink()
cat("\nDone. Report and figures are in", normalizePath(out_dir), "\n")
