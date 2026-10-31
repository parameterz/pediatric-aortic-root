# =====================================================================
#  Aortic root z-score models: SOV, STJ and AAo
#  Same pipeline as the AAo work, run for each structure:
#    exponents -> height-only cubic -> location-scale GAM
#    (height, weight-for-height, age, sex) -> z checks -> export tables
#
#  Writes:
#    aao_explore/root_model.json      tables for the web calculator
#    aao_explore/root_testcases.csv   R z-scores to check the JavaScript
#    aao_explore/09_local_exponents.png
#  Needs: mgcv, jsonlite, ggplot2
# =====================================================================

library(mgcv)
library(jsonlite)
library(ggplot2)

out_dir <- "aao_explore"
dir.create(out_dir, showWarnings = FALSE)
raw <- echo_zscores_phn

# ---- 0. SET THESE: column names for each structure (values in cm) ----
print(grep("^AR_", names(raw), value = TRUE))     # find the right names here
# Use whichever spelling exists (with or without the underscore)
pick <- function(...) { o <- intersect(c(...), names(raw)); if (length(o)) o[1] else c(...)[1] }
structures <- c(SOV = pick("AR_ARD", "ARARD"),     # aortic root, sinuses of Valsalva
                STJ = pick("AR_ASJD", "ARASJD"),   # sinotubular junction
                AAO = pick("AR_AAD", "ARAAD"))     # ascending aorta
labels <- c(SOV = "Sinuses of Valsalva", STJ = "Sinotubular junction", AAO = "Ascending aorta")
stopifnot("Set the structure column names above" = all(structures %in% names(raw)))

# Rows judged implausible (height/weight) in the AAo work
drop_prob <- c("1957", "1908", "2627", "647", "2136")
# Structure-specific measurement errors (excluded for that structure only)
drop_extra <- list(STJ = c("1080"),        # 17.9 y, 172 cm, STJ 13 mm (z -6.6)
                   AAO = c("489", "236"))  # ~140 cm, AAo 12 and 13 mm (z -4.7, -4.0); normal SOV and STJ

# ---- 1. Shared body-size data and weight-for-height reference --------
base <- data.frame(age = raw$ECHO_AGE, ht = raw$ECHO_HT / 100, wt = raw$ECHO_WT,
                   male = as.integer(grepl("^m", tolower(as.character(raw$SEX)))))
for (k in names(structures)) base[[k]] <- raw[[structures[[k]]]]
base <- base[complete.cases(base[, c("age", "ht", "wt", "male")]), ]
base$lht  <- log(base$ht)
base$sexf <- factor(base$male, 0:1, c("F", "M"))
base$bmi  <- base$wt / base$ht^2

wfh_fit <- gam(log(wt) ~ s(lht, k = 10), data = base, method = "REML")
base$wfh <- log(base$wt) - predict(wfh_fit, base)

cat("\nUnits check (should look like cm):\n")
print(sapply(names(structures), function(k) summary(base[[k]])))

skew <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
zsum <- function(z, d) {
  sl <- function(v) { x <- d[[v]]; if (v != "male") x <- (x - mean(x)) / sd(x); unname(coef(lm(z ~ x))[2]) }
  dec  <- tapply(z, cut(d$ht, quantile(d$ht, 0:10 / 10), include.lowest = TRUE), mean)
  decs <- tapply(z, cut(d$ht, quantile(d$ht, 0:10 / 10), include.lowest = TRUE), sd)
  agem <- tapply(z, cut(d$age, c(-Inf, 1/12, 3, 6, 12, 16, Inf), right = FALSE), mean)
  round(c(sd = sd(z), skew = skew(z), pct_gt2 = 100 * mean(z > 2), pct_lt_m2 = 100 * mean(z < -2),
          slope_lht = sl("lht"), slope_age = sl("age"), slope_wt = sl("wt"), slope_bmi = sl("bmi"),
          slope_male = sl("male"),
          decile_mean_min = min(dec), decile_mean_max = max(dec),
          decile_sd_min = min(decs), decile_sd_max = max(decs),
          agegrp_mean_min = min(agem), agegrp_mean_max = max(agem)), 3)
}

# Shared grids (tables clamp outside the data range)
lht_grid <- seq(min(base$lht), max(base$lht), length.out = 401)
wfh_rng  <- range(base$wfh[!rownames(base) %in% drop_prob])
wfh_grid <- seq(wfh_rng[1], wfh_rng[2], length.out = 201)
age_grid <- seq(0, 18, length.out = 181)
elw      <- unname(predict(wfh_fit, data.frame(lht = lht_grid)))
interp   <- function(x, xs, ys) approx(xs, ys, x, rule = 2)$y

# ---- 2. Fit, check and tabulate one structure -------------------------
fit_one <- function(key) {
  d <- base[!is.na(base[[key]]) & base[[key]] > 0 &
            !rownames(base) %in% c(drop_prob, drop_extra[[key]]), ]
  d$lad <- log(d[[key]])
  cat("\n\n================ ", labels[[key]], " (", key, "), n = ", nrow(d), " ================\n", sep = "")

  ex <- function(s) { m <- lm(lad ~ lht, s); ci <- confint(m)[2, ]
                      c(n = nrow(s), exponent = coef(m)[[2]], lo95 = ci[[1]], hi95 = ci[[2]]) }
  cat("\nHeight exponents (log diameter ~ log height):\n")
  print(round(rbind(all = ex(d), `<2 y` = ex(d[d$age < 2, ]), `>=2 y` = ex(d[d$age >= 2, ])), 3))

  # Height + sex cubic (fallback when weight is missing; sex is always known)
  cub <- lm(lad ~ poly(lht, 3, raw = TRUE) + male, d)
  d$z_cub <- resid(cub) / sigma(cub)

  # Full location-scale GAM
  g <- gam(list(lad ~ s(lht, k = 10) + s(wfh, k = 6) + s(age, k = 6) + sexf, ~ s(lht, k = 6)),
           family = gaulss(), data = d, method = "REML", select = TRUE)
  cat("\nGAM smooth terms:\n"); print(round(summary(g)$s.table, 3))
  cat("Male coefficient (log scale):", round(coef(g)[["sexfM"]], 4), "\n")

  pred <- function(nd) {
    p <- predict(g, nd, type = "response")
    list(mu = unname(p[, 1]), sd = unname(if (median(p[, 2]) > 1) 1 / p[, 2] else p[, 2]))
  }
  pf <- pred(d)
  d$z <- (d$lad - pf$mu) / pf$sd

  cat("\nz checks (ideal: sd 1, skew 0, ~2.3% each tail, slopes 0, decile/age means near 0):\n")
  print(rbind(gam = zsum(d$z, d), cubic = zsum(d$z_cub, d)))

  cat("\nLargest |z| (check for measurement or entry errors):\n")
  o <- order(-abs(d$z))[1:8]
  print(transform(d[o, c("age", "ht", "wt", key, "z")], ht = ht * 100, z = round(z, 2)))

  # Tables
  lev <- c("F", "M")
  ref <- data.frame(lht = 0, wfh = 0, age = 0, sexf = factor("F", levels = lev))
  grid_df <- function(var, vals) { x <- ref[rep(1, length(vals)), ]; x[[var]] <- vals; x }
  mu_ref <- pred(ref)$mu
  p_lht <- pred(grid_df("lht", lht_grid))
  tab <- list(
    label = labels[[key]], n = nrow(d),
    f_lht = p_lht$mu, sd = p_lht$sd,
    f_wfh = pred(grid_df("wfh", wfh_grid))$mu - mu_ref,
    f_age = pred(grid_df("age", age_grid))$mu - mu_ref,
    male  = pred(transform(ref, sexf = factor("M", levels = lev)))$mu - mu_ref,
    cubic = list(b = unname(coef(cub)[1:4]), male = unname(coef(cub)[5]), sd = sigma(cub))
  )

  # Validate tables against the model
  zt <- (d$lad - (interp(d$lht, lht_grid, tab$f_lht) + interp(d$wfh, wfh_grid, tab$f_wfh) +
                  interp(d$age, age_grid, tab$f_age) + tab$male * d$male)) /
        interp(d$lht, lht_grid, tab$sd)
  cat("\nMax |z(table) - z(model)|:", signif(max(abs(zt - d$z)), 3), "\n")

  # Test cases (z computed exactly as the JavaScript will)
  set.seed(match(key, names(structures)))
  s <- d[sample(nrow(d), 10), ]
  tc <- data.frame(structure = key, diam_mm = s[[key]] * 10, height_cm = s$ht * 100,
                   weight_kg = s$wt, age_y = s$age, sex = as.character(s$sexf),
                   z_R = round(zt[match(rownames(s), rownames(d))], 4))
  list(tab = tab, tc = tc, d = d, cub = cub)
}

fits <- lapply(names(structures), fit_one)
names(fits) <- names(structures)

# ---- 3. Local exponent by structure (from the height-only cubic) -------
le <- do.call(rbind, lapply(names(fits), function(k) {
  b <- coef(fits[[k]]$cub); L <- lht_grid
  data.frame(structure = labels[[k]], ht_cm = exp(L) * 100, exponent = b[2] + 2 * b[3] * L + 3 * b[4] * L^2)
}))
p <- ggplot(le, aes(ht_cm, exponent, colour = structure)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") + geom_hline(yintercept = 0.5, linetype = "dotted") +
  labs(x = "Height (cm)", y = "Local exponent: d log(diameter) / d log(height)", colour = NULL,
       title = "Local height exponent by structure", subtitle = "From the height-only cubic; dashed = 1, dotted = 0.5")
ggsave(file.path(out_dir, "09_local_exponents.png"), p, width = 8, height = 5, dpi = 150); print(p)

# ---- 4. Chaudhary root (SOV) check: supplement vs refit ----------------
if ("SOV" %in% names(fits)) {
  dS <- fits$SOV$d
  dS$csa_ht <- pi * (dS$SOV / 2)^2 / dS$ht
  cat("\n\nCSA/ht aortic root: refit on these data vs the supplement spreadsheet\n")
  print(rbind(refit = round(coef(lm(csa_ht ~ age * male, dS)), 4),
              sheet = c(1.6414, 0.1294, 0.0962, 0.0147)))
}

# ---- 5. Export ----------------------------------------------------------
model_json <- list(
  meta = list(
    name = "Aortic root z-scores (location-scale GAMs)",
    source = "NIH/NHLBI Pediatric Heart Network Echo-Z/Normal ECG public use dataset, downloaded 01/01/2026",
    created = as.character(Sys.Date()),
    form = "z = (ln(D_cm) - mu) / sd; mu = f_lht(L) + f_wfh(wfh) + f_age(age) + male*isMale; sd = sd(L); L = ln(height_m); wfh = ln(weight_kg) - elw(L)",
    ranges = list(height_cm = round(100 * exp(range(lht_grid)), 1), wfh = round(wfh_rng, 3), age = c(0, 18))
  ),
  lht = lht_grid, elw = elw, wfh = wfh_grid, age = age_grid,
  structures = lapply(fits, `[[`, "tab")
)
write_json(model_json, file.path(out_dir, "root_model.json"), auto_unbox = TRUE, digits = 8)
tcs <- do.call(rbind, lapply(fits, `[[`, "tc"))
write.csv(tcs, file.path(out_dir, "root_testcases.csv"), row.names = FALSE)
cat("\nWrote root_model.json (", file.size(file.path(out_dir, "root_model.json")), " bytes) and ",
    nrow(tcs), " test cases.\n", sep = "")
