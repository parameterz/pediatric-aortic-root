# =====================================================================
#  Distribution diagnostics for the BCCG models (SOV, STJ, AAo)
#
#  Run after aortic_root_models.R in the same session (uses `fits`, `labels`, `out_dir`).
#    1. Box-Cox power (lambda): free vs fixed at 0 (log), 1 (linear), 2 (squared),
#       and BCT (adds a kurtosis parameter) for comparison
#    2. Royston-Wright Q-tests by height group: mean, variance, skewness, kurtosis
#    3. Worm plots by height group, saved as aao_explore/11_worm_<level>.png
#    4. Moments and tail rates of the z-scores
#  Loops run at top level on purpose: gamlss looks up data by name.
# =====================================================================

library(gamlss)

ctrl   <- gamlss.control(n.cyc = 200, trace = FALSE)
skew   <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
exkurt <- function(x) { x <- x - mean(x); mean(x^4) / mean(x^2)^2 - 3 }
pct    <- function(z, c) c(100 * mean(z > c), 100 * mean(z < -c))

lam_tab <- list(); diag_tab <- list()

for (k in c("SOV", "STJ", "AAO")) {
  dk <- fits[[k]]$d[, c(k, "lht", "wfh", "age", "sexf")]
  names(dk)[1] <- "y"
  f_mu <- y ~ pb(lht) + pb(wfh) + pb(age) + sexf

  m <- gamlss(f_mu, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
              family = BCCG(mu.link = "log"), data = dk, control = ctrl)

  # ---- 1. Shape: free lambda vs fixed scales, and BCT ----
  m_log <- update(m, nu.fix = TRUE, nu.start = 0.001)
  m_lin <- update(m, nu.fix = TRUE, nu.start = 1)
  m_sq  <- update(m, nu.fix = TRUE, nu.start = 2)
  m_bct <- try(gamlss(f_mu, sigma.formula = ~ pb(lht), nu.formula = ~ 1, tau.formula = ~ 1,
                      family = BCT(mu.link = "log"), data = dk, control = ctrl, start.from = m),
               silent = TRUE)
  se <- suppressWarnings(vcov(m, type = "se"))
  lam_tab[[k]] <- c(lambda = unname(coef(m, "nu")), se = unname(se[length(se)]),
                    AIC_free = AIC(m), AIC_log = AIC(m_log), AIC_linear = AIC(m_lin),
                    AIC_squared = AIC(m_sq),
                    AIC_BCT = if (inherits(m_bct, "try-error")) NA else AIC(m_bct))

  # ---- 2. Q-tests by height group ----
  cat("\n==== ", labels[[k]], ": Q-tests by height group ====\n", sep = "")
  print(Q.stats(m, xvar = dk$lht, n.inter = 8, plot = FALSE))
  cat("Note: the Z1 (mean) Q-test subtracts all mu degrees of freedom, including the\n",
      "weight-for-height and age smooths, so its p-value is conservative.\n", sep = "")

  # ---- 3. Worm plots ----
  png(file.path(out_dir, paste0("11_worm_", k, ".png")), width = 1600, height = 1100, res = 150)
  wp(m, xvar = dk$lht, n.inter = 6, ylim.worm = 0.6)
  dev.off()

  # ---- 4. Moments and tails ----
  z  <- resid(m)
  sw <- shapiro.test(z)
  diag_tab[[k]] <- c(mean = mean(z), sd = sd(z), skew = skew(z), ex_kurt = exkurt(z),
                     setNames(pct(z, 2),   c("gt2", "lt_m2")),
                     setNames(pct(z, 2.5), c("gt2.5", "lt_m2.5")),
                     setNames(pct(z, 3),   c("gt3", "lt_m3")),
                     shapiro_W = unname(sw$statistic), shapiro_p = sw$p.value)
}

cat("\nShape (lambda) and family comparison (lower AIC = better):\n")
print(round(do.call(rbind, lam_tab), 3))
cat("\nMoments and tail rates of z:\n")
print(round(do.call(rbind, diag_tab), 3))
cat("Expected % beyond each side: 2 -> 2.28, 2.5 -> 0.62, 3 -> 0.13\n")
