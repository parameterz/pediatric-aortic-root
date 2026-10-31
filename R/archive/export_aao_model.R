# =====================================================================
#  Export the final AAo GAM (b_prob) as lookup tables for a web calculator
#
#  Run after phn_aao_explore.R and the multivariable section, so these
#  exist: d, wfh_fit, b_prob, drop_prob
#
#  Writes:
#    aao_explore/aao_model.json   -- tables the JavaScript reads
#    aao_explore/aao_testcases.csv -- R z-scores to check the JS against
# =====================================================================

library(mgcv)
library(jsonlite)      # install.packages("jsonlite") if needed

m     <- b_prob
d_fit <- d[!rownames(d) %in% drop_prob, ]
lev   <- c("F", "M")

pred <- function(nd) {
  p  <- predict(m, nd, type = "response")
  sd <- if (median(p[, 2]) > 1) 1 / p[, 2] else p[, 2]   # col 2 may be 1/sd
  list(mu = unname(p[, 1]), sd = unname(sd))
}

# ---- Grids over the observed data range --------------------------------
lht_grid <- seq(min(d_fit$lht), max(d_fit$lht), length.out = 401)
wfh_grid <- seq(min(d_fit$wfh), max(d_fit$wfh), length.out = 201)
age_grid <- seq(0, 18, length.out = 181)

# Reference child: lht = 0 (1 m), wfh = 0, age = 0, female
ref <- data.frame(lht = 0, wfh = 0, age = 0, sexf = factor("F", levels = lev))
grid_df <- function(var, vals) { g <- ref[rep(1, length(vals)), ]; g[[var]] <- vals; g }
mu_ref <- pred(ref)$mu

# Additive model, so each term can be tabulated separately:
#   f_lht includes the intercept (other terms at their reference values)
#   f_wfh, f_age are 0 at the reference point
p_lht <- pred(grid_df("lht", lht_grid))
f_lht <- p_lht$mu
sd_lht <- p_lht$sd                                   # SD depends on lht only
f_wfh <- pred(grid_df("wfh", wfh_grid))$mu - mu_ref
f_age <- pred(grid_df("age", age_grid))$mu - mu_ref
male  <- pred(transform(ref, sexf = factor("M", levels = lev)))$mu - mu_ref

# Weight-for-height reference: expected log weight for a given log height
elw <- unname(predict(wfh_fit, data.frame(lht = lht_grid)))

# ---- Validate: table + linear interpolation vs the real model ----------
interp <- function(x, xs, ys) approx(xs, ys, x, rule = 2)$y
tab_z <- function(aad_cm, ht_m, wt_kg, age, is_male) {
  L   <- log(ht_m)
  wfh <- log(wt_kg) - interp(L, lht_grid, elw)
  mu  <- interp(L, lht_grid, f_lht) + interp(wfh, wfh_grid, f_wfh) +
         interp(age, age_grid, f_age) + male * is_male
  (log(aad_cm) - mu) / interp(L, lht_grid, sd_lht)
}
p_fit   <- pred(d_fit)
z_model <- (d_fit$lad - p_fit$mu) / p_fit$sd
z_table <- tab_z(d_fit$aad, d_fit$ht, d_fit$wt, d_fit$age, d_fit$male)
cat("Max |z(table) - z(model)| over", nrow(d_fit), "kids:",
    signif(max(abs(z_table - z_model)), 3), "\n")
cat("Max |wfh(table) - wfh(fit)|:",
    signif(max(abs(log(d_fit$wt) - interp(d_fit$lht, lht_grid, elw) - d_fit$wfh)), 3), "\n")

# ---- Export --------------------------------------------------------------
model_json <- list(
  meta = list(
    name    = "PHN ascending aorta z-score (location-scale GAM)",
    source  = "PHN Normal Echocardiogram Z-Score Database, healthy non-obese children <= 18 y",
    n       = nrow(d_fit),
    created = as.character(Sys.Date()),
    form    = "z = (ln(AAD_cm) - mu) / sd;  mu = f_lht(L) + f_wfh(wfh) + f_age(age) + male*isMale;  sd = sd(L);  L = ln(height_m);  wfh = ln(weight_kg) - elw(L)",
    ranges  = list(height_cm = round(100 * exp(range(d_fit$lht)), 1),
                   wfh = round(range(d_fit$wfh), 3),
                   age = c(0, 18))
  ),
  lht = lht_grid, f_lht = f_lht, sd = sd_lht, elw = elw,
  wfh = wfh_grid, f_wfh = f_wfh,
  age = age_grid, f_age = f_age,
  male = male
)
write_json(model_json, file.path(out_dir, "aao_model.json"),
           auto_unbox = TRUE, digits = 8, pretty = FALSE)

# Test cases for checking the JavaScript, including Dan's 8 y girl
set.seed(1)
tc <- d_fit[sample(nrow(d_fit), 25),
            c("aad", "ht", "wt", "age", "male")]
tc <- rbind(tc, data.frame(aad = 2.4, ht = 1.28, wt = 28, age = 8, male = 0,
                           row.names = "dan_case"))
tc$z_R <- tab_z(tc$aad, tc$ht, tc$wt, tc$age, tc$male)
tc <- data.frame(aad_mm = tc$aad * 10, height_cm = tc$ht * 100, weight_kg = tc$wt,
                 age_y = tc$age, sex = ifelse(tc$male == 1, "M", "F"),
                 z_R = round(tc$z_R, 4), row.names = rownames(tc))
write.csv(tc, file.path(out_dir, "aao_testcases.csv"))
print(tc["dan_case", ])
cat("JSON size:", file.size(file.path(out_dir, "aao_model.json")), "bytes\n")
