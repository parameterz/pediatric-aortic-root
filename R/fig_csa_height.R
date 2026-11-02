# =====================================================================
#  Figure: is CSA/height independent of height?
#
#  Run after aortic_root_models.R in the same session (uses fits, out_dir).
#  Chaudhary et al. (Pediatr Cardiol 2026) plotted CSA/height against age. A size-
#  neutral index should be flat against body size, so the plot that tests the index
#  is CSA/height against HEIGHT. This script draws it for the aortic root (sinuses)
#  and the ascending aorta, with:
#    - a smooth (GAM) of CSA/height on height, with 95% band
#    - the power-law exponent b from ln(CSA/ht) = a + b ln(ht), with 95% CI
#      (b = 0 for a size-neutral index; diameter exponent = (b + 1) / 2)
#    - the share of variance in CSA/height explained by height alone (R^2)
#
#  Two versions:
#    fig_csa_height_scatter.(pdf|png|tiff)   every child as a point (like the paper's Figure 1)
#    fig_csa_height_binned.(pdf|png)         20 equal-count height groups: mean and middle 50%
#                                            (aggregated only; suitable for the web site)
#  Check the PHN data-use terms before publishing the scatter version.
# =====================================================================

library(ggplot2)
library(mgcv)

levels_fig <- c(SOV = "Aortic root (sinuses of Valsalva)", AAO = "Ascending aorta")
okabe <- c(Female = "#D55E00", Male = "#0072B2")

d_all <- do.call(rbind, lapply(names(levels_fig), function(k) {
  d <- fits[[k]]$d
  data.frame(level = levels_fig[[k]], ht_cm = d$ht * 100, ht = d$ht, age = d$age,
             sex = ifelse(d$male == 1, "Male", "Female"),
             csa_ht = pi * (d[[k]] / 2)^2 / d$ht)                     # cm^2 / m, as in the paper
}))
d_all$level <- factor(d_all$level, levels = levels_fig)

# ---- Statistics per level ----------------------------------------------------
stats <- do.call(rbind, lapply(levels(d_all$level), function(lv) {
  d  <- d_all[d_all$level == lv, ]
  pw <- lm(log(csa_ht) ~ log(ht), d)
  ci <- confint(pw)[2, ]
  r2_ht  <- summary(lm(csa_ht ~ poly(ht_cm, 3), d))$r.squared            # height alone (flexible)
  r2_age <- summary(lm(csa_ht ~ age * sex, d))$r.squared                  # the paper's model form
  r2_both <- summary(lm(csa_ht ~ poly(ht_cm, 3) + age * sex, d))$r.squared
  data.frame(level = lv, n = nrow(d), b = coef(pw)[2], lo = ci[1], hi = ci[2],
             diam_exp = (coef(pw)[2] + 1) / 2, r2_ht = r2_ht, r2_age = r2_age, r2_both = r2_both)
}))
stats$level <- factor(stats$level, levels = levels_fig)
cat("\n==== CSA/height against height ====\n")
print(transform(stats, b = round(b, 3), lo = round(lo, 3), hi = round(hi, 3), diam_exp = round(diam_exp, 3),
                r2_ht = round(r2_ht, 3), r2_age = round(r2_age, 3), r2_both = round(r2_both, 3)), row.names = FALSE)
cat("b = 0 would mean CSA/height is independent of height (the premise of the index).\n",
    "r2_age: variance explained by the paper's model (age x sex); r2_ht: by height alone; r2_both: height + age x sex.\n",
    "If r2_ht >= r2_age and adding age adds little, the 'age effect' is mostly residual height dependence.\n", sep = "")

stats$label <- sprintf("CSA/ht ~ height^%.2f (95%% CI %.2f to %.2f)\nimplied diameter exponent %.2f\nheight explains %.0f%% of variance",
                       stats$b, stats$lo, stats$hi, stats$diam_exp, 100 * stats$r2_ht)

# smooth for each level
grid <- do.call(rbind, lapply(levels(d_all$level), function(lv) {
  d <- d_all[d_all$level == lv, ]
  g <- gam(csa_ht ~ s(ht_cm, k = 10), data = d, method = "REML")
  nd <- data.frame(ht_cm = seq(min(d$ht_cm), max(d$ht_cm), length.out = 200))
  p <- predict(g, nd, se.fit = TRUE)
  data.frame(level = lv, ht_cm = nd$ht_cm, fit = p$fit, lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit)
}))
grid$level <- factor(grid$level, levels = levels_fig)

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11), legend.position = "bottom")
ymax <- max(d_all$csa_ht) * 1.02

# ---- Version 1: scatter ----------------------------------------------------------
p1 <- ggplot(d_all, aes(ht_cm, csa_ht)) +
  geom_point(aes(colour = sex), size = 0.7, alpha = 0.35, stroke = 0) +
  geom_ribbon(data = grid, aes(ht_cm, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey20", alpha = 0.2) +
  geom_line(data = grid, aes(ht_cm, fit), inherit.aes = FALSE, linewidth = 0.9, colour = "black") +
  geom_text(data = stats, aes(x = 45, y = ymax, label = label), hjust = 0, vjust = 1, size = 3, lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~ level) +
  scale_colour_manual(values = okabe, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = "Height (cm)", y = expression("Cross-sectional area / height (cm"^2*"/m)")) +
  base_theme

# ---- Version 2: binned (aggregated only) -------------------------------------------
bins <- do.call(rbind, lapply(levels(d_all$level), function(lv) {
  d  <- d_all[d_all$level == lv, ]
  g  <- cut(d$ht_cm, unique(quantile(d$ht_cm, 0:20 / 20)), include.lowest = TRUE)
  data.frame(level = lv, x = tapply(d$ht_cm, g, median), m = tapply(d$csa_ht, g, mean),
             q25 = tapply(d$csa_ht, g, quantile, 0.25), q75 = tapply(d$csa_ht, g, quantile, 0.75),
             n = as.vector(table(g)))
}))
bins$level <- factor(bins$level, levels = levels_fig)
p2 <- ggplot(bins, aes(x, m)) +
  geom_ribbon(data = grid, aes(ht_cm, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey20", alpha = 0.15) +
  geom_line(data = grid, aes(ht_cm, fit), inherit.aes = FALSE, linewidth = 0.8, colour = "grey30") +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0, colour = "#0072B2", linewidth = 0.7) +
  geom_point(colour = "#0072B2", size = 2) +
  geom_text(data = stats, aes(x = 45, y = ymax * 0.8, label = label), hjust = 0, vjust = 1, size = 3, lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~ level) +
  coord_cartesian(ylim = c(0, ymax * 0.8)) +
  labs(x = "Height (cm)", y = expression("Cross-sectional area / height (cm"^2*"/m)"),
       caption = "Points: mean of 20 equal-count height groups (about 160 children each); bars: middle 50%;\nline: smooth with 95% band.") +
  base_theme + theme(plot.caption = element_text(hjust = 0, size = 8), plot.caption.position = "plot")

ggsave(file.path(out_dir, "fig_csa_height_scatter.pdf"),  p1, width = 7, height = 3.8)
ggsave(file.path(out_dir, "fig_csa_height_scatter.png"),  p1, width = 7, height = 3.8, dpi = 300)
ggsave(file.path(out_dir, "fig_csa_height_scatter.tiff"), p1, width = 7, height = 3.8, dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "fig_csa_height_binned.pdf"),   p2, width = 7, height = 3.8)
ggsave(file.path(out_dir, "fig_csa_height_binned.png"),   p2, width = 7, height = 3.8, dpi = 300)
cat("\nWrote fig_csa_height_scatter (pdf/png/tiff) and fig_csa_height_binned (pdf/png) to", out_dir, "\n")
print(p1)
