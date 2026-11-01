# =====================================================================
#  Head-to-head: every normalization method applied to the same healthy PHN children
#
#  Run after aortic_root_models.R in the same session (uses base, fits, labels, out_dir).
#  Needs the gamlss package and the Lauzon-Schnittka authors' zscorecalculator.r (ls_file).
#
#  If a method is a good reference for these children, its z-scores should have
#  mean 0 and SD 1, about 2.3% beyond each of +/-2, no skew, and no trend left
#  against height, age, build or sex. Methods compared, per level:
#    ours_cv    this model (BCCG), 10-fold cross-validated: each child scored by a model
#               fitted without them, so no home advantage
#    ours       this model, in sample (what the calculator gives)
#    hs         the height-and-sex model (no weight), in sample
#    phn        PHN Echo Z-Score Project table, D / BSA^0.5 (Haycock), all ages
#    mah_bsa    Mahgerefteh 2021, D / BSA^0.5 (>= 2 y; SOV and AAo)
#    mah_ht     Mahgerefteh 2021, D / height (>= 2 y; SOV and AAo)
#    ls         Lauzon-Schnittka 2025, the authors' own z_score() (clinical cohort; not refit)
#    csa_tab    Chaudhary 2026 CSA/ht, age-group mean and SD (SOV and AAo)
#    csa_sheet  Chaudhary 2026 CSA/ht, continuous age, coefficients as distributed
#    csa_refit  Chaudhary 2026 CSA/ht, continuous age, coefficients refit to these data
#  (The PHN, Mahgerefteh and Chaudhary methods were built on these same PHN children,
#   so they also have a home advantage; Lauzon-Schnittka is the only external one.)
#
#  Exports AGGREGATED summaries only (no individual records):
#    - moments, tails and trend slopes per level and method
#    - mean and SD of z in 20 equal-count bins of height, age, weight for height, BMIz
#    - agreement on z > 2 with the cross-validated model (counts)
#    - expected diameter vs BMI-for-age z for three example children (model-to-model)
#  Writes aao_explore/methods.json and methods.js (window.AORTIC_METHODS); copy both to models/.
#  Loops that call gamlss run at top level on purpose (gamlss looks up data by name).
# =====================================================================

library(gamlss)
library(jsonlite)

ls_file <- "zscorecalculator.r"          # path to the Lauzon-Schnittka authors' script
k_folds <- 10
set.seed(20260101)

ctrl   <- gamlss.control(n.cyc = 200, trace = FALSE)
lev    <- c("F", "M")
bccg_z <- function(y, mu, sigma, nu) ifelse(abs(nu) < 1e-6, log(y / mu) / sigma, ((y / mu)^nu - 1) / (nu * sigma))
haycock <- function(ht_cm, wt) 0.024265 * ht_cm^0.3964 * wt^0.5378

# ---- Published constants (same as js/aortic-model.js) ------------------------
PHN <- list(SOV = c(2.06, 0.18), STJ = c(1.69, 0.16), AAO = c(1.79, 0.18))             # cm / m (BSA^0.5)
MAH <- list(SOV = list(ht = c(1.61, 0.15), bsa = c(2.07, 0.18)),
            AAO = list(ht = c(1.39, 0.14), bsa = c(1.78, 0.18)))
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

# ---- Lauzon-Schnittka: the authors' function, plus their WHO BMIz -------------
source(ls_file)
ls_env <- new.env()
for (s in as.list(body(z_score))[-1])
  if (is.call(s) && identical(s[[1]], as.name("<-")) && is.name(s[[2]]) &&
      (grepl("^bmiz_[lms]_[mf]", as.character(s[[2]])) || as.character(s[[2]]) %in% c("interpolate", "bmiz_who")))
    eval(s, ls_env)
who_bmiz <- function(ht_cm, wt, age, sx) mapply(ls_env$bmiz_who, ht_cm, wt, age, sx)
ls_z     <- function(lv, d_mm, sx, ht_cm, wt, age) z_score(lv, d_mm, sx, ht_cm, wt, age)

# ---- 1. Score every child under every method ----------------------------------
scores <- list()
cvfit  <- list()     # full-data models for section 4
for (k in names(fits)) {
  cat("\n==== ", labels[[k]], " ====\n", sep = "")
  d  <- fits[[k]]$d
  dk <- d[, c(k, "lht", "wfh", "age", "sexf")]; names(dk)[1] <- "y"
  ht_cm <- d$ht * 100; d_mm <- d[[k]] * 10; sx <- ifelse(d$male == 1, "m", "f")

  # this model and the height-and-sex model, in sample
  m_all <- gamlss(y ~ pb(lht) + pb(wfh) + pb(age) + sexf, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
                  family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  m_hs  <- gamlss(y ~ pb(lht) + sexf, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
                  family = BCCG(mu.link = "log"), data = dk, control = ctrl)
  cvfit[[k]] <- list(m = m_all, dk = dk)
  z_ours <- resid(m_all); z_hs <- resid(m_hs)

  # 10-fold cross-validation (folds stratified by height so each fold spans all sizes)
  fold <- integer(nrow(dk))
  fold[order(dk$lht)] <- rep_len(sample(k_folds), nrow(dk))
  z_cv <- rep(NA_real_, nrow(dk))
  for (f in seq_len(k_folds)) {
    tr <- dk[fold != f, ]
    te <- dk[fold == f, ]
    mf <- gamlss(y ~ pb(lht) + pb(wfh) + pb(age) + sexf, sigma.formula = ~ pb(lht), nu.formula = ~ 1,
                 family = BCCG(mu.link = "log"), data = tr, control = ctrl)
    # keep held-out predictors inside the training range (smooths are not extrapolated)
    te_c <- te
    for (v in c("lht", "wfh", "age")) te_c[[v]] <- pmin(pmax(te[[v]], min(tr[[v]])), max(tr[[v]]))
    pa <- predictAll(mf, newdata = te_c[, c("lht", "wfh", "age", "sexf")], data = tr)
    z_cv[fold == f] <- bccg_z(te$y, pa$mu, pa$sigma, pa$nu)
    cat("  fold", f, "converged:", mf$converged, "\n")
  }

  # published indices
  bsa_rt <- sqrt(haycock(ht_cm, d$wt))
  z_phn  <- (d[[k]] / bsa_rt - PHN[[k]][1]) / PHN[[k]][2]
  ge2    <- d$age >= 2
  z_mah_bsa <- z_mah_ht <- rep(NA_real_, nrow(d))
  if (!is.null(MAH[[k]])) {
    z_mah_bsa[ge2] <- ((d[[k]] / bsa_rt - MAH[[k]]$bsa[1]) / MAH[[k]]$bsa[2])[ge2]
    z_mah_ht[ge2]  <- ((d[[k]] / d$ht   - MAH[[k]]$ht[1])  / MAH[[k]]$ht[2])[ge2]
  }

  # Lauzon-Schnittka, the authors' code (slow: one child at a time)
  cat("  Lauzon-Schnittka z_score() on", nrow(d), "children ...\n")
  z_ls <- ls_z(k, d_mm, sx, ht_cm, d$wt, d$age)

  # Chaudhary CSA / height
  z_csa_tab <- z_csa_sheet <- z_csa_refit <- rep(NA_real_, nrow(d))
  if (!is.null(CSA[[k]])) {
    cs   <- CSA[[k]]
    csa  <- pi * (d[[k]] / 2)^2 / d$ht                           # cm^2 / m
    bin  <- findInterval(d$age, cs$lo)
    tabm <- ifelse(d$male == 1, cs$M[bin, 1], cs$F[bin, 1])
    tabs <- ifelse(d$male == 1, cs$M[bin, 2], cs$F[bin, 2])
    reg  <- function(b) b[["int"]] + b[["age"]] * d$age + d$male * (b[["male"]] + b[["ma"]] * d$age)
    z_csa_tab   <- (csa - tabm) / tabs
    z_csa_sheet <- (csa - reg(cs$sheet)) / tabs                  # SD from the age-group table, as in the calculator
    z_csa_refit <- (csa - reg(cs$refit)) / tabs
  }

  scores[[k]] <- data.frame(
    ht_cm = ht_cm, age = d$age, wfh = 100 * (exp(d$wfh) - 1), male = d$male,
    bmiz = who_bmiz(ht_cm, d$wt, d$age, sx),
    ours_cv = z_cv, ours = z_ours, hs = z_hs, phn = z_phn, mah_bsa = z_mah_bsa, mah_ht = z_mah_ht,
    ls = z_ls, csa_tab = z_csa_tab, csa_sheet = z_csa_sheet, csa_refit = z_csa_refit)
}

methods <- c(ours_cv = "This model, cross-validated", ours = "This model, in sample", hs = "Height and sex",
             phn = "PHN z-score (BSA^0.5)", mah_bsa = "Mahgerefteh BSA^0.5 (>= 2 y)", mah_ht = "Mahgerefteh height (>= 2 y)",
             ls = "Lauzon-Schnittka", csa_tab = "CSA/ht, age-group table",
             csa_sheet = "CSA/ht, continuous age (as distributed)", csa_refit = "CSA/ht, continuous age (refit)")

# ---- 2. Summaries ---------------------------------------------------------------
skew <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
kurt <- function(x) { x <- x - mean(x); mean(x^4) / mean(x^2)^2 - 3 }
slope_sd <- function(z, x) unname(coef(lm(z ~ I((x - mean(x)) / sd(x))))[2])   # change in z per SD of x

summ <- list(); bins <- list(); agree <- list()
bin_vars <- c(ht_cm = "Height (cm)", age = "Age (y)", wfh = "Weight for height (%)", bmiz = "BMI-for-age z (WHO)")
for (k in names(scores)) {
  S <- scores[[k]]
  ms <- names(methods)[sapply(names(methods), function(m) sum(!is.na(S[[m]])) > 50)]
  summ[[k]] <- do.call(rbind, lapply(ms, function(m) {
    ok <- !is.na(S[[m]]) & is.finite(S[[m]]); z <- S[[m]][ok]; s <- S[ok, ]
    data.frame(method = m, n = sum(ok), mean = mean(z), sd = sd(z), skew = skew(z), exkurt = kurt(z),
               gt2 = 100 * mean(z > 2), lt2 = 100 * mean(z < -2), gt3 = 100 * mean(z > 3), lt3 = 100 * mean(z < -3),
               sl_ht = slope_sd(z, s$ht_cm), sl_age = slope_sd(z, s$age), sl_wfh = slope_sd(z, s$wfh),
               sl_bmiz = slope_sd(z, s$bmiz), sex_diff = mean(z[s$male == 1]) - mean(z[s$male == 0]))
  }))
  cat("\n==== ", labels[[k]], ": z-score summary by method (target: mean 0, SD 1, 2.3% beyond each of +/-2, slopes 0) ====\n", sep = "")
  print(format(transform(summ[[k]], method = methods[method]), digits = 2, nsmall = 2), row.names = FALSE)

  # binned trends: 20 equal-count bins of each body-size variable
  bins[[k]] <- lapply(names(bin_vars), function(v) {
    x  <- S[[v]]
    br <- unique(quantile(x, 0:20 / 20, na.rm = TRUE))
    g  <- cut(x, br, include.lowest = TRUE)
    out <- list(var = v, label = bin_vars[[v]],
                x = unname(tapply(x, g, median)), xlo = unname(tapply(x, g, min)), xhi = unname(tapply(x, g, max)),
                n = unname(as.vector(table(g))))
    out$mean <- lapply(setNames(ms, ms), function(m) unname(tapply(S[[m]], g, function(z) if (sum(!is.na(z)) >= 20) mean(z, na.rm = TRUE) else NA)))
    out$sd   <- lapply(setNames(ms, ms), function(m) unname(tapply(S[[m]], g, function(z) if (sum(!is.na(z)) >= 20) sd(z, na.rm = TRUE) else NA)))
    out
  })

  # agreement with the cross-validated model on z > 2
  agree[[k]] <- do.call(rbind, lapply(setdiff(ms, "ours_cv"), function(m) {
    ok <- !is.na(S[[m]]) & !is.na(S$ours_cv)
    a <- S$ours_cv[ok] > 2; b <- S[[m]][ok] > 2
    data.frame(method = m, n = sum(ok), both = sum(a & b), ours_only = sum(a & !b), other_only = sum(!a & b), neither = sum(!a & !b))
  }))
  cat("\nAgreement on z > 2 with the cross-validated model (counts):\n")
  print(transform(agree[[k]], method = methods[method]), row.names = FALSE)
}

# ---- 3. Model to model: expected diameter vs BMI-for-age z ----------------------
# Three example children at a fixed height and age; weight set from the WHO BMI-for-age
# z (-2 to +4). Expected (median) diameter under each weight-aware method, as % of its
# value at BMIz 0. Our model holds weight for height flat beyond the PHN range, flagged.
examples <- data.frame(id = c("boy3", "boy10", "girl14"), label = c("Boy, 3 y, 95 cm", "Boy, 10 y, 138 cm", "Girl, 14 y, 160 cm"),
                       sex = c("m", "m", "f"), age = c(3, 10, 14), ht_cm = c(95, 138, 160))
bz_grid <- seq(-2, 4, by = 0.25)
wt_for_bmiz <- function(bz, ht_cm, age, sx)
  uniroot(function(w) ls_env$bmiz_who(ht_cm, w, age, sx) - bz, c(3, 250), tol = 1e-8)$root
ls_median <- function(lv, sx, ht_cm, wt, age)
  uniroot(function(dm) z_score(lv, dm, sx, ht_cm, wt, age), c(2, 80), tol = 1e-8)$root

curves <- list()
for (k in names(cvfit)) {
  m_all <- cvfit[[k]]$m
  dk    <- cvfit[[k]]$dk
  wrng  <- range(dk$wfh)
  curves[[k]] <- lapply(seq_len(nrow(examples)), function(i) {
    e  <- examples[i, ]
    wt <- sapply(bz_grid, wt_for_bmiz, ht_cm = e$ht_cm, age = e$age, sx = e$sex)
    wfh_raw <- log(wt) - unname(predict(wfh_fit, data.frame(lht = rep(log(e$ht_cm / 100), length(wt)))))
    nd <- data.frame(lht = log(e$ht_cm / 100), wfh = pmin(pmax(wfh_raw, wrng[1]), wrng[2]), age = e$age,
                     sexf = factor(ifelse(e$sex == "m", "M", "F"), levels = lev))
    ours <- predictAll(m_all, newdata = nd, data = dk)$mu * 10                     # mm
    lsm  <- sapply(wt, function(w) ls_median(k, e$sex, e$ht_cm, w, e$age))
    phn  <- PHN[[k]][1] * sqrt(haycock(e$ht_cm, wt)) * 10
    rel  <- function(v) 100 * v / approx(bz_grid, v, 0)$y
    list(id = e$id, label = e$label, bmiz = bz_grid, wt = wt, wfh_pct = 100 * (exp(wfh_raw) - 1),
         beyond_range = wfh_raw < wrng[1] | wfh_raw > wrng[2],
         mm = list(ours = ours, ls = lsm, phn = phn), pct = list(ours = rel(ours), ls = rel(lsm), phn = rel(phn)))
  })
  cat("\n", labels[[k]], ": expected diameter at BMIz +3 as % of BMIz 0 (ours / L-S / PHN BSA):\n", sep = "")
  for (cc in curves[[k]]) {
    j <- which.min(abs(cc$bmiz - 3))
    cat(sprintf("  %-20s %5.1f%% / %5.1f%% / %5.1f%%%s\n", cc$label, cc$pct$ours[j], cc$pct$ls[j], cc$pct$phn[j],
                if (cc$beyond_range[j]) "   (ours held flat: beyond PHN build range)" else ""))
  }
}

# ---- 4. Export --------------------------------------------------------------------
out <- list(meta = list(created = as.character(Sys.Date()), k_folds = k_folds,
                        note = "Aggregated comparison of normalization methods on the PHN Echo-Z cohort; no individual records"),
            methods = as.list(methods), labels = as.list(labels),
            summary = summ, bins = bins, agree = agree, curves = curves)
txt <- toJSON(out, auto_unbox = TRUE, digits = 5, na = "null")
writeLines(txt, file.path(out_dir, "methods.json"))
writeLines(paste0("window.AORTIC_METHODS = ", txt, ";"), file.path(out_dir, "methods.js"))
cat("\nWrote methods.json / methods.js. Copy both to models/.\n")
