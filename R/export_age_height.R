# =====================================================================
#  Height-for-age reference from the PHN cohort, for the calculator's age check
#
#  Run after aortic_root_models.R in the same session (uses base, drop_prob, out_dir).
#  Fits height ~ age separately for girls and boys with BCCG (LMS-type):
#    median and spread as smooth functions of sqrt(age), constant shape.
#  Exports the median, sigma and lambda on an age grid (0 to 18 y by 0.05 y):
#  summary curves only, no individual records.
#
#  The web pages use it to flag a height that is far from typical for the entered
#  age (|z| > 4), for example an age entered in months instead of years.
#
#  Writes aao_explore/age_height.json and aao_explore/age_height.js; copy both to models/.
#  Loops run at top level on purpose (gamlss looks up data by name).
# =====================================================================

library(gamlss)
library(jsonlite)

ctrl     <- gamlss.control(n.cyc = 200, trace = FALSE)
age_grid <- seq(0, 18, by = 0.05)

dh <- base[!rownames(base) %in% drop_prob, c("ht", "age", "sexf")]
dh$ht <- dh$ht * 100           # cm
dh$sa <- sqrt(dh$age)          # sqrt(age) straightens the steep infant growth

out <- list(age = age_grid)
for (s in c("F", "M")) {
  ds <- dh[dh$sexf == s, c("ht", "sa")]
  m  <- gamlss(ht ~ pb(sa), sigma.formula = ~ pb(sa), nu.formula = ~ 1,
               family = BCCG, data = ds, control = ctrl)
  p  <- predictAll(m, newdata = data.frame(sa = sqrt(age_grid)), data = ds)
  z  <- resid(m)
  cat(sprintf("%s: n = %d, lambda = %.2f, |z| > 4: %d children, max |z| = %.2f\n",
              s, nrow(ds), coef(m, "nu"), sum(abs(z) > 4), max(abs(z))))
  out[[s]] <- list(mu = unname(p$mu), sigma = unname(p$sigma), nu = unname(coef(m, "nu")))
}

# Quick look: median height at a few ages
chk <- sapply(c("F", "M"), function(s) round(approx(age_grid, out[[s]]$mu, c(0, 1, 2, 5, 8, 12, 15, 18))$y, 1))
rownames(chk) <- paste0(c(0, 1, 2, 5, 8, 12, 15, 18), " y")
print(chk)

txt <- toJSON(c(list(meta = list(created = as.character(Sys.Date()),
                                 note = "Height-for-age reference (BCCG on sqrt(age)) from the PHN Echo-Z cohort; summary curves only")),
                out), auto_unbox = TRUE, digits = 6)
writeLines(txt, file.path(out_dir, "age_height.json"))
writeLines(paste0("window.AORTIC_AGEHT = ", txt, ";"), file.path(out_dir, "age_height.js"))
cat("Wrote age_height.json / age_height.js. Copy both to models/.\n")
