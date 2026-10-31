# =====================================================================
#  Export model diagnostics for the web page (aggregated summaries only)
#
#  Run after aortic_root_models.R in the same session (uses fits, labels, out_dir).
#  For each level, refits the published BCCG model and summarizes its z-scores:
#    - z by bins of height, age, weight, BMI and weight for height (20 bins of ~160 children:
#      mean, SD, quantiles, tail rates) plus a smoothed trend with its standard error
#    - z by sex
#    - Royston-Wright Q-tests by height group
#    - Q-Q points at fixed probabilities, worm-plot deviations by height group, histogram counts
#    - moments, tail rates, Shapiro-Wilk
#    - shape comparison (lambda free vs fixed at 0, 1, 2; BCT)
#    - local height exponent of the diameter (smooth derivative with SE)
#  No individual records are written: every number summarizes a group of children.
#
#  Writes aao_explore/diagnostics.json and aao_explore/diagnostics.js
#  Copy both to models/. Loops run at top level on purpose (gamlss looks up data by name).
# =====================================================================

library(gamlss)
library(mgcv)
library(jsonlite)

ctrl   <- gamlss.control(n.cyc = 200, trace = FALSE)
skew   <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
exkurt <- function(x) { x <- x - mean(x); mean(x^4) / mean(x^2)^2 - 3 }
probs_qq   <- c(0.005, 0.01, 0.02, seq(0.05, 0.95, by = 0.025), 0.98, 0.99, 0.995)
probs_worm <- seq(0.02, 0.98, by = 0.02)

binstats <- function(x, z, nb = 20) {
  br <- unique(quantile(x, seq(0, 1, length.out = nb + 1)))
  g  <- cut(x, br, include.lowest = TRUE)
  out <- do.call(rbind, lapply(split(seq_along(z), g), function(i) {
    zz <- z[i]
    q  <- quantile(zz, c(.025, .1, .25, .5, .75, .9, .975), names = FALSE)
    data.frame(x = median(x[i]), xlo = min(x[i]), xhi = max(x[i]), n = length(i),
               mean = mean(zz), sd = sd(zz),
               q025 = q[1], q10 = q[2], q25 = q[3], q50 = q[4], q75 = q[5], q90 = q[6], q975 = q[7],
               gt2 = 100 * mean(zz > 2), lt2 = 100 * mean(zz < -2))
  }))
  rownames(out) <- NULL
  out
}
trend <- function(x, z, n = 60) {
  g  <- gam(z ~ s(x, k = 10), method = "REML")
  xs <- seq(quantile(x, .01), quantile(x, .99), length.out = n)
  p  <- predict(g, data.frame(x = xs), se.fit = TRUE)
  data.frame(x = xs, fit = unname(p$fit), se = unname(p$se.fit))
}

levels_out <- list()

for (k in c("SOV", "STJ", "AAO")) {
  cat("\n==== ", labels[[k]], " ====\n", sep = "")
  d0 <- fits[[k]]$d
  dk <- d0[, c(k, "lht", "wfh", "age", "sexf")]
  names(dk)[1] <- "y"
  f_mu <- y ~ pb(lht) + pb(wfh) + pb(age) + sexf

  m <- gamlss(f_mu, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
              family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  z <- resid(m)

  # ---- shape comparison ----
  m_log <- update(m, nu.fix = TRUE, nu.start = 0.001)
  m_lin <- update(m, nu.fix = TRUE, nu.start = 1)
  m_sq  <- update(m, nu.fix = TRUE, nu.start = 2)
  m_bct <- try(gamlss(f_mu, sigma.formula = ~ pb(lht), nu.formula = ~ 1, tau.formula = ~ 1,
                      family = BCT(mu.link = "log"), data = dk, control = ctrl, start.from = m),
               silent = TRUE)
  se <- suppressWarnings(vcov(m, type = "se"))
  aic <- list(free = AIC(m), log = AIC(m_log), linear = AIC(m_lin), squared = AIC(m_sq),
              bct = if (inherits(m_bct, "try-error")) NA else AIC(m_bct))

  # ---- covariates ----
  cov <- list(height = exp(d0$lht) * 100, age = d0$age, weight = d0$wt, bmi = d0$bmi,
              wfh = 100 * (exp(d0$wfh) - 1))
  bins   <- lapply(cov, binstats, z = z)
  trends <- lapply(cov, trend, z = z)
  sex <- do.call(rbind, lapply(split(z, d0$sexf), function(zz)
    data.frame(n = length(zz), mean = mean(zz), sd = sd(zz), se = sd(zz) / sqrt(length(zz)),
               gt2 = 100 * mean(zz > 2), lt2 = 100 * mean(zz < -2))))
  sex$sex <- rownames(sex); rownames(sex) <- NULL

  # ---- Q-tests ----
  qs <- Q.stats(m, xvar = dk$lht, n.inter = 8, plot = FALSE)
  qs <- data.frame(group = rownames(qs), as.data.frame(qs), check.names = FALSE, row.names = NULL)
  # label height groups in cm instead of log height
  lab_cm <- function(s) {
    v <- suppressWarnings(as.numeric(strsplit(trimws(s), "\\s+to\\s+")[[1]]))
    if (length(v) == 2 && all(is.finite(v))) sprintf("%.0f–%.0f cm", 100 * exp(v[1]), 100 * exp(v[2])) else s
  }
  qs$group <- vapply(qs$group, lab_cm, "")

  # ---- distribution ----
  qq <- data.frame(p = probs_qq, theo = qnorm(probs_qq), samp = unname(quantile(z, probs_qq)))
  hb <- seq(floor(min(z) * 4) / 4, ceiling(max(z) * 4) / 4, by = 0.25)
  h  <- hist(z, breaks = hb, plot = FALSE)
  grp <- cut(dk$lht, unique(quantile(dk$lht, seq(0, 1, length.out = 7))), include.lowest = TRUE)
  worm <- lapply(split(seq_along(z), grp), function(i) {
    hh <- exp(range(dk$lht[i])) * 100
    list(label = sprintf("%.0f–%.0f cm", hh[1], hh[2]), n = length(i),
         theo = qnorm(probs_worm), dev = unname(quantile(z[i], probs_worm)) - qnorm(probs_worm))
  })
  names(worm) <- NULL

  # ---- local exponent (height-only smooth on the log scale) ----
  g  <- gam(lad ~ s(lht, k = 10), data = d0, method = "REML")
  L  <- seq(min(d0$lht), max(d0$lht), length.out = 120); eps <- 1e-4
  Xd <- (predict(g, data.frame(lht = L + eps), type = "lpmatrix") -
         predict(g, data.frame(lht = L - eps), type = "lpmatrix")) / (2 * eps)
  expo <- data.frame(ht = exp(L) * 100, b = drop(Xd %*% coef(g)),
                     se = sqrt(rowSums((Xd %*% vcov(g)) * Xd)))

  sw <- shapiro.test(z)
  levels_out[[k]] <- list(
    label = labels[[k]], n = length(z),
    lambda = unname(coef(m, "nu")), lambda_se = unname(se[length(se)]),
    moments = list(mean = mean(z), sd = sd(z), skew = skew(z), ex_kurt = exkurt(z),
                   shapiro_W = unname(sw$statistic), shapiro_p = sw$p.value),
    tails = data.frame(cut = c(2, 2.5, 3), expected = 100 * pnorm(-c(2, 2.5, 3)),
                       above = sapply(c(2, 2.5, 3), function(c) 100 * mean(z > c)),
                       below = sapply(c(2, 2.5, 3), function(c) 100 * mean(z < -c))),
    aic = aic, bins = bins, trends = trends, sex = sex, qstats = qs,
    qq = qq, hist = list(breaks = h$breaks, counts = h$counts), worm = worm, exponent = expo)
  cat("n =", length(z), "| lambda =", round(coef(m, "nu"), 3), "| done\n")
}

diag <- list(
  meta = list(created = as.character(Sys.Date()),
              source = "NIH/NHLBI Pediatric Heart Network Echo-Z/Normal ECG public use dataset, downloaded 01/01/2026",
              note = "Aggregated summaries only; no individual records."),
  levels = levels_out)
txt <- toJSON(diag, dataframe = "columns", auto_unbox = TRUE, digits = 5, na = "null")
writeLines(txt, file.path(out_dir, "diagnostics.json"))
writeLines(paste0("window.AORTIC_DIAG = ", txt, ";"), file.path(out_dir, "diagnostics.js"))
cat("\nWrote diagnostics.json / diagnostics.js (", file.size(file.path(out_dir, "diagnostics.json")),
    " bytes). Copy both to models/.\n", sep = "")
