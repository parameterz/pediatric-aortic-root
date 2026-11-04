# =====================================================================
#  Does an age-based reference stand in for height? Height-for-age bias
#
#  Run after aortic_root_models.R in the same session (uses base, fits, drop_prob, out_dir).
#  Needs gamlss.
#
#  Idea: CSA/height still rises with height. Modeling it against age corrects that only
#  for a child of average height for age. A child tall for age (e.g., Marfan syndrome)
#  should then score high, and a child short for age (e.g., Turner syndrome) low,
#  even with an aorta normal for their size.
#
#  Steps
#   1. Height-for-age z (HAZ) for every child, from an LMS-type (BCCG) reference fitted to
#      the PHN children themselves, by sex (same approach as export_age_height.R).
#   2. z-scores under each method (as in export_methods_comparison.R):
#        csa_tab    Chaudhary CSA/ht, age-group mean and SD
#        csa_sheet  Chaudhary CSA/ht, continuous age, coefficients as distributed
#        csa_refit  Chaudhary CSA/ht, continuous age, coefficients refit
#        phn        PHN z-score, D / BSA^0.5 (size-based, no age)
#        mah_ht     Mahgerefteh D / height (>= 2 y; size-based, no age)
#        ours       this model (height, weight for height, age, sex), in sample
#   3. Within-age effect: lm(z ~ HAZ) (HAZ is already age- and sex-standardized), with the
#      predicted shift in z for a child at +2 SD (Marfan-like) and -2.5 SD (Turner-like),
#      plus mean z and % above +2 / below -2 in HAZ groups.
#  Output: console tables and an aggregated plot (HAZ groups only; no individual records).
#  Loops that call gamlss run at top level on purpose (gamlss looks up data by name).
# =====================================================================

library(gamlss)
library(ggplot2)

ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)

# ---- 1. Height-for-age z from the PHN children ---------------------------------------
hb <- base[!rownames(base) %in% drop_prob, c("ht", "age", "sexf")]
hb$ht_cm <- hb$ht * 100
hb$sa    <- sqrt(hb$age)
hb$haz   <- NA_real_
for (s in c("F", "M")) {
  hs <- hb[hb$sexf == s, c("ht_cm", "sa")]
  mh <- gamlss(ht_cm ~ pb(sa), sigma.formula = ~ pb(sa), nu.formula = ~ 1,
               family = BCCG, data = hs, control = ctrl)
  hb$haz[hb$sexf == s] <- resid(mh)                 # normalized quantile residuals = z
  cat(sprintf("Height-for-age reference (%s): n = %d, SD of HAZ = %.2f\n", s, nrow(hs), sd(resid(mh))))
}

# ---- 2. z under each method ---------------------------------------------------------
haycock <- function(ht_cm, wt) 0.024265 * ht_cm^0.3964 * wt^0.5378
PHN <- list(SOV = c(2.06, 0.18), AAO = c(1.79, 0.18))
MAH <- list(SOV = c(1.61, 0.15), AAO = c(1.39, 0.14))
CSA <- list(
  SOV = list(lo = c(0, 1/12, 3, 6, 12, 16),
             F = rbind(c(1.4128447, 0.2550527), c(1.8518291, 0.4176144), c(2.2326539, 0.3698699),
                       c(2.5486409, 0.4543347), c(2.9936960, 0.4995857), c(3.1964701, 0.5267062)),
             M = rbind(c(1.5346930, 0.2775356), c(1.9816039, 0.4120444), c(2.4048436, 0.3998567),
                       c(2.8695985, 0.5095964), c(3.2347010, 0.5950158), c(3.6397332, 0.6259371)),
             sheet = c(int = 1.6414159903, age = 0.1294167470, male = 0.0961634070, ma = 0.0147000314),
             refit = c(int = 1.6413, age = 0.0962, male = 0.1302, ma = 0.0146)),
  AAO = list(lo = c(0, 1/12, 3, 6, 12, 16),
             F = rbind(c(1.1566680, 0.2718176), c(1.3894141, 0.3164411), c(1.6674859, 0.3513082),
                       c(1.9498509, 0.3534838), c(2.3370879, 0.4802758), c(2.4963250, 0.5387956)),
             M = rbind(c(1.2500483, 0.2862416), c(1.3956703, 0.3201912), c(1.7176683, 0.3666948),
                       c(2.0984815, 0.4530524), c(2.3480720, 0.4751761), c(2.6417850, 0.5027709)),
             sheet = c(int = 1.2612332551, age = 0.0538404716, male = 0.0751494789, ma = 0.0031328790),
             refit = c(int = 1.2612, age = 0.0751, male = 0.0553, ma = 0.0031)))

methods <- c(csa_tab = "CSA/ht, age-group table", csa_sheet = "CSA/ht, continuous (as distributed)",
             csa_refit = "CSA/ht, continuous (refit)", phn = "PHN z-score (BSA^0.5)",
             mah_ht = "Mahgerefteh height (>= 2 y)", ours = "This model (in sample)")

zz <- list()
for (k in c("SOV", "AAO")) {
  d  <- fits[[k]]$d
  dk <- d[, c(k, "lht", "wfh", "age", "sexf")]; names(dk)[1] <- "y"
  m_ours <- gamlss(y ~ pb(lht) + pb(wfh) + pb(age) + sexf, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
                   family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  cs   <- CSA[[k]]
  csa  <- pi * (d[[k]] / 2)^2 / d$ht
  bin  <- findInterval(d$age, cs$lo)
  tabm <- ifelse(d$male == 1, cs$M[bin, 1], cs$F[bin, 1])
  tabs <- ifelse(d$male == 1, cs$M[bin, 2], cs$F[bin, 2])
  reg  <- function(b) b[["int"]] + b[["age"]] * d$age + d$male * (b[["male"]] + b[["ma"]] * d$age)
  ge2  <- d$age >= 2
  zz[[k]] <- data.frame(
    haz = hb[rownames(d), "haz"], age = d$age,
    csa_tab = (csa - tabm) / tabs, csa_sheet = (csa - reg(cs$sheet)) / tabs, csa_refit = (csa - reg(cs$refit)) / tabs,
    phn = (d[[k]] / sqrt(haycock(d$ht * 100, d$wt)) - PHN[[k]][1]) / PHN[[k]][2],
    mah_ht = ifelse(ge2, (d[[k]] / d$ht - MAH[[k]][1]) / MAH[[k]][2], NA),
    ours = resid(m_ours))
}

# ---- 3. Within-age effect of height for age -----------------------------------------
haz_grp <- function(h) cut(h, c(-Inf, -2, -1, 1, 2, Inf), labels = c("< -2", "-2 to -1", "-1 to +1", "+1 to +2", "> +2"))
res <- list(); grp <- list()
for (k in names(zz)) {
  Z <- zz[[k]]
  res[[k]] <- do.call(rbind, lapply(names(methods), function(m) {
    ok <- !is.na(Z[[m]]) & !is.na(Z$haz)
    f  <- lm(Z[[m]][ok] ~ Z$haz[ok]); ci <- confint(f)[2, ]
    data.frame(level = k, method = methods[[m]], n = sum(ok), slope = coef(f)[2], lo = ci[1], hi = ci[2],
               shift_tall = 2 * coef(f)[2], shift_short = -2.5 * coef(f)[2])
  }))
  g <- haz_grp(Z$haz)
  grp[[k]] <- do.call(rbind, lapply(names(methods), function(m) {
    do.call(rbind, lapply(levels(g), function(L) {
      z <- Z[[m]][g == L & !is.na(Z[[m]])]
      data.frame(level = k, method = methods[[m]], haz_group = L, n = length(z),
                 mean_z = if (length(z) >= 10) mean(z) else NA,
                 above2 = if (length(z) >= 10) 100 * mean(z > 2) else NA,
                 below2 = if (length(z) >= 10) 100 * mean(z < -2) else NA)
    }))
  }))
}
res <- do.call(rbind, res); grp <- do.call(rbind, grp)
rownames(res) <- NULL

cat("\n==== Shift in z per SD of height for age (within age) ====\n")
cat("shift_tall: predicted z shift at +2 SD height for age (Marfan-like); shift_short: at -2.5 SD (Turner-like)\n")
print(transform(res, slope = round(slope, 3), lo = round(lo, 3), hi = round(hi, 3),
                shift_tall = round(shift_tall, 2), shift_short = round(shift_short, 2)), row.names = FALSE)
cat("\n==== Mean z and tails by height-for-age group ====\n")
print(transform(grp, mean_z = round(mean_z, 2), above2 = round(above2, 1), below2 = round(below2, 1)), row.names = FALSE)

# ---- Plot (aggregated: HAZ groups) --------------------------------------------------
grp$level <- factor(ifelse(grp$level == "SOV", "Aortic root (sinuses)", "Ascending aorta"), levels = c("Aortic root (sinuses)", "Ascending aorta"))
grp$method <- factor(grp$method, levels = methods)
grp$haz_group <- factor(grp$haz_group, levels = levels(haz_grp(0)))
p <- ggplot(grp[!is.na(grp$mean_z), ], aes(haz_group, mean_z, colour = method, group = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~ level) +
  scale_colour_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#009E73", "#0072B2"), name = NULL) +
  labs(x = "Height for age (z)", y = "Mean z-score of aortic diameter",
       caption = "Healthy PHN children grouped by height for age. A reference that adjusts for size stays at 0 across groups.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot") +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "fig_height_for_age_bias.png"), p, width = 8, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "fig_height_for_age_bias.pdf"), p, width = 8, height = 4.5)
cat("\nWrote fig_height_for_age_bias (png/pdf) to", out_dir, "\n")
print(p)
