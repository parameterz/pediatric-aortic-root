# =====================================================================
#  Export BCCG (LMS-type) models for SOV, STJ and AAo as lookup tables
#
#  Run after aortic_root_models.R in the same session; uses:
#    fits, lht_grid, wfh_grid, age_grid, elw, wfh_rng, labels
#
#  Model per structure (BCCG = Box-Cox Cole-Green):
#    median  log(mu) = f(log height) + g(weight for height) + h(age) + sex   (log link)
#    spread  sigma = s(log height)                                        (log link)
#    shape   nu    = constant  (Box-Cox power, lambda)
#    z = ((y/mu)^nu - 1) / (nu * sigma)        [log(y/mu)/sigma when nu = 0]
#  Fallback without weight: log(mu) = f(log height) + sex, same sigma/nu structure.
#
#  Writes aao_explore/root_bccg.json and aao_explore/root_bccg_testcases.csv
#  Note: loops run at top level on purpose (gamlss looks up data by name).
# =====================================================================

library(gamlss)
library(jsonlite)

ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)
lev  <- c("F", "M")
mk   <- function(lht = 0, wfh = 0, age = 0, sex = "F")
  data.frame(lht = lht, wfh = wfh, age = age, sexf = factor(sex, levels = lev))
interp <- function(x, xs, ys) approx(xs, ys, x, rule = 2)$y
bccg_z <- function(y, mu, sigma, nu)
  if (abs(nu) < 1e-6) log(y / mu) / sigma else ((y / mu)^nu - 1) / (nu * sigma)

out <- list(); tcs <- list()

for (k in names(fits)) {
  cat("\n==== ", labels[[k]], " ====\n", sep = "")
  dk <- fits[[k]]$d[, c(k, "lht", "wfh", "age", "sexf")]
  names(dk)[1] <- "y"

  m  <- gamlss(y ~ pb(lht) + pb(wfh) + pb(age) + sexf, sigma.formula = ~ pb(lht),
               nu.formula = ~ 1, family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  m0 <- gamlss(y ~ pb(lht) + sexf, sigma.formula = ~ pb(lht),
               nu.formula = ~ 1, family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  cat("Converged:", m$converged, m0$converged,
      "| lambda full =", round(coef(m, "nu"), 3), "| lambda fallback =", round(coef(m0, "nu"), 3), "\n")

  # ---- tables (full model; log(mu) is additive) ----
  pr <- predictAll(m, newdata = mk(),                data = dk)
  pl <- predictAll(m, newdata = mk(lht = lht_grid),  data = dk)
  pw <- predictAll(m, newdata = mk(wfh = wfh_grid),  data = dk)
  pa <- predictAll(m, newdata = mk(age = age_grid),  data = dk)
  pm <- predictAll(m, newdata = mk(sex = "M"),       data = dk)
  full <- list(nu = unname(coef(m, "nu")),
               mu_lht = unname(log(pl$mu)), sigma = unname(pl$sigma),
               mu_wfh = unname(log(pw$mu) - log(pr$mu)), mu_age = unname(log(pa$mu) - log(pr$mu)),
               male = unname(log(pm$mu) - log(pr$mu)))

  # ---- tables (height + sex fallback) ----
  qr <- predictAll(m0, newdata = mk(),               data = dk)
  ql <- predictAll(m0, newdata = mk(lht = lht_grid), data = dk)
  qm <- predictAll(m0, newdata = mk(sex = "M"),      data = dk)
  hs <- list(nu = unname(coef(m0, "nu")), mu_lht = unname(log(ql$mu)), sigma = unname(ql$sigma),
             male = unname(log(qm$mu) - log(qr$mu)))

  # ---- validate tables against gamlss's own z (normalized quantile residuals) ----
  male <- as.integer(dk$sexf == "M")
  mu_t <- exp(interp(dk$lht, lht_grid, full$mu_lht) + interp(dk$wfh, wfh_grid, full$mu_wfh) +
              interp(dk$age, age_grid, full$mu_age) + full$male * male)
  z_t  <- bccg_z(dk$y, mu_t, interp(dk$lht, lht_grid, full$sigma), full$nu)
  mu_h <- exp(interp(dk$lht, lht_grid, hs$mu_lht) + hs$male * male)
  z_h  <- bccg_z(dk$y, mu_h, interp(dk$lht, lht_grid, hs$sigma), hs$nu)
  cat("Max |mu(table) - mu(model)| (cm):", signif(max(abs(mu_t - fitted(m, "mu"))), 3), "\n")
  cat("Max |z(table) - z(gamlss)|  full:", signif(max(abs(z_t - resid(m))), 3),
      " fallback:", signif(max(abs(z_h - resid(m0))), 3), "\n")
  cat("Tails (full): >+2 =", round(100 * mean(z_t > 2), 2), "%, <-2 =", round(100 * mean(z_t < -2), 2), "%\n")

  out[[k]] <- list(label = labels[[k]], n = nrow(dk), full = full, hs = hs)

  set.seed(match(k, names(fits)) + 10)
  i <- sample(nrow(dk), 10)
  tcs[[k]] <- data.frame(structure = k, diam_mm = dk$y[i] * 10, height_cm = round(exp(dk$lht[i]) * 100, 2),
                         weight_kg = fits[[k]]$d$wt[i], age_y = dk$age[i], sex = as.character(dk$sexf[i]),
                         z_R = round(z_t[i], 4), z_R_hs = round(z_h[i], 4))
}

model_json <- list(
  meta = list(
    name = "Aortic root z-scores (BCCG / LMS-type GAMLSS)",
    source = "NIH/NHLBI Pediatric Heart Network Echo-Z/Normal ECG public use dataset, downloaded 01/01/2026",
    created = as.character(Sys.Date()),
    form = "z = ((D/mu)^nu - 1)/(nu*sigma); ln mu (cm) = mu_lht(L) + mu_wfh(wfh) + mu_age(age) + male*isMale; sigma = sigma(L); L = ln(height_m); wfh = ln(weight_kg) - elw(L)",
    ranges = list(height_cm = round(100 * exp(range(lht_grid)), 1), wfh = round(wfh_rng, 3), age = c(0, 18))
  ),
  lht = lht_grid, elw = elw, wfh = wfh_grid, age = age_grid,
  structures = out
)
write_json(model_json, file.path(out_dir, "root_bccg.json"), auto_unbox = TRUE, digits = 8)
tc <- do.call(rbind, tcs)
write.csv(tc, file.path(out_dir, "root_bccg_testcases.csv"), row.names = FALSE)
cat("\nWrote root_bccg.json (", file.size(file.path(out_dir, "root_bccg.json")), " bytes) and ",
    nrow(tc), " test cases.\n", sep = "")

# JavaScript wrappers used by the web pages (classic scripts, so pages also open straight from disk)
writeLines(paste0("window.AORTIC_MODEL = ", toJSON(model_json, auto_unbox = TRUE, digits = 8), ";"),
           file.path(out_dir, "root_bccg.js"))
writeLines(paste0("window.AORTIC_TESTS = ", toJSON(tc, dataframe = "rows", digits = 8), ";"),
           file.path(out_dir, "root_bccg_testcases.js"))
cat("Copy root_bccg.json and root_bccg.js to models/, and root_bccg_testcases.csv/.js to tests/.\n")
