# =====================================================================
#  Lauzon-Schnittka 2025 model for the web pages, taken from the authors' own R code
#
#  Lauzon-Schnittka J, et al. Circ Cardiovasc Imaging. 2025;18:e017944.
#  The authors distribute an R function, z_score(), in zscorecalculator.r. This script
#  sources that file, reads the coefficients, spline knots and weights and the WHO
#  BMI-for-age LMS vectors straight out of the function body (nothing is retyped),
#  and writes them for the browser. It then runs z_score() itself on a set of test
#  children so the browser engine can be checked against the authors' code.
#
#  Set ls_file to your copy of zscorecalculator.r (not included in this repository).
#  Stand-alone: does not need the PHN data.
#
#  Writes (to out_dir; copy to the repo as noted):
#    lauzon_schnittka.json / .js    -> models/   (window.LS_MODEL)
#    ls_reference_cases.js          -> tests/    (window.LS_TESTS)
# =====================================================================

library(jsonlite)

ls_file <- "zscorecalculator.r"           # path to the authors' script
if (!exists("out_dir")) out_dir <- "aao_explore"
dir.create(out_dir, showWarnings = FALSE)

source(ls_file)                           # defines z_score()
stopifnot(exists("z_score"))

# ---- 1. Pull the objects defined inside z_score() --------------------------
# Every top-level assignment in the function body whose name we need is evaluated
# in a scratch environment: the WHO vectors, the parameter list and its order.
keep <- function(nm) grepl("^bmiz_[lms]_[mf]", nm) || nm %in% c("parameters", "echoordervar", "interpolate", "bmiz_who")
env  <- new.env()
for (s in as.list(body(z_score))[-1]) {
  if (is.call(s) && identical(s[[1]], as.name("<-")) && is.name(s[[2]]) && keep(as.character(s[[2]])))
    eval(s, env)
}

levels <- list()
for (lv in c("SOV", "STJ", "AAO")) for (sx in c("m", "f")) {
  key <- paste0(lv, sx)
  el  <- env$parameters[[match(key, env$echoordervar)]]
  spl <- function(i) if (is.null(el[[i]]) || is.null(el[[i + 1]])) NULL else list(t = el[[i]], w = el[[i + 1]])
  levels[[key]] <- list(p = el[[1]], b = unlist(el[2:10]),
                        s1 = spl(11), s2 = spl(13), s3 = spl(15), s4 = spl(17))
  cat(sprintf("%s: p = %.2f, splines: %s\n", key, el[[1]],
              paste(c("s1 mu~age", "s2 mu~BMIz", "s3 sigma~age", "s4 lambda~age")[!sapply(c(11, 13, 15, 17), function(i) is.null(el[[i]]))], collapse = ", ")))
}

who <- function(sx) list(
  day = list(L = env[[paste0("bmiz_l_", sx, "_2y")]], M = env[[paste0("bmiz_m_", sx, "_2y")]], S = env[[paste0("bmiz_s_", sx, "_2y")]]),
  mon = list(L = env[[paste0("bmiz_l_", sx)]],        M = env[[paste0("bmiz_m_", sx)]],        S = env[[paste0("bmiz_s_", sx)]]))
bmiz_tab <- list(m = who("m"), f = who("f"))
max_age  <- (length(bmiz_tab$m$mon$L) + 23) / 12   # the monthly table interpolates up to (not including) this age
cat(sprintf("WHO tables: %d daily rows (< 2 y), %d monthly rows (2 y to %.2f y)\n",
            length(bmiz_tab$m$day$L), length(bmiz_tab$m$mon$L), max_age))

model <- list(levels = levels, bmiz = bmiz_tab, max_age = max_age,
              meta = list(created = as.character(Sys.Date()),
                          source = "Lauzon-Schnittka J, et al. Z scores for pediatric echocardiography dimensions adjusted for body size, BMI, and age. Circ Cardiovasc Imaging. 2025;18:e017944",
                          note = "Coefficients, B-spline knots and weights, and WHO BMI-for-age LMS vectors read from the authors' R function z_score() (zscorecalculator.r)"))
txt <- toJSON(model, auto_unbox = TRUE, digits = NA, null = "null")
writeLines(txt, file.path(out_dir, "lauzon_schnittka.json"))
writeLines(paste0("window.LS_MODEL = ", txt, ";"), file.path(out_dir, "lauzon_schnittka.js"))

# ---- 2. Reference cases from z_score() itself ------------------------------
bmiz1 <- function(ht, wt, age, sx) env$bmiz_who(ht, wt, age, sx)
d_at  <- function(lv, sx, ht, wt, age, z) {     # diameter at a given z, by inverting z_score()
  f <- function(d) z_score(lv, d, sx, ht, wt, age) - z
  uniroot(f, c(0.5, 120), tol = 1e-12)$root
}

set.seed(2025)
ages <- c(0.001, 0.01, 0.05, 0.25, 0.5, 1, 1.5, 1.99, 2, 2.01, 3, 5, 7, 9, 11, 13, 15, 16.5, 18, 18.9)
cases <- do.call(rbind, lapply(seq_along(ages), function(i) {
  a  <- ages[i]
  ht <- round(50 + 25 * a^0.6 * (1 + 0.03 * sin(i)) + 6 * a^1.1 / (1 + a^0.4), 1)    # rough height for age
  bmi_f <- c(0.8, 1, 1.25, 1.6, 2.3)[(i %% 5) + 1]                                    # lean to very obese (BMIz beyond the clamp)
  bmi0  <- if (a < 1) 14 + 3 * a else 16 + 0.25 * max(0, a - 6)
  wt <- round(bmi0 * bmi_f * (ht / 100)^2, 2)
  do.call(rbind, lapply(c("m", "f"), function(sx) {
    k  <- 1 + 0.12 * sin(3 * i + (sx == "f"))        # spread the diameters around typical
    data.frame(sex = sx, age_y = a, ht = ht, wt = wt,
               SOV = round(k * 0.9 * ht^0.75, 1), STJ = round(k * 0.75 * ht^0.75, 1), AAO = round(k * 0.8 * ht^0.75 * 1.05, 1))
  }))
}))
# Very obese and very lean extras, plus an exact-integer age
cases <- rbind(cases,
  data.frame(sex = c("m", "f", "f", "m"), age_y = c(10, 10, 4, 12), ht = c(140, 140, 102, 150), wt = c(70, 70, 12.5, 30),
             SOV = c(27, 25, 17, 26), STJ = c(22, 21, 14, 21), AAO = c(26, 25, 16, 24)))

tests <- lapply(seq_len(nrow(cases)), function(i) {
  r <- cases[i, ]
  z  <- sapply(c("SOV", "STJ", "AAO"), function(lv) z_score(lv, r[[lv]], r$sex, r$ht, r$wt, r$age_y))
  d0 <- sapply(c("SOV", "STJ", "AAO"), function(lv) d_at(lv, r$sex, r$ht, r$wt, r$age_y, 0))
  d2 <- sapply(c("SOV", "STJ", "AAO"), function(lv) d_at(lv, r$sex, r$ht, r$wt, r$age_y, 2))
  list(sex = r$sex, age_y = r$age_y, ht = r$ht, wt = r$wt, SOV = r$SOV, STJ = r$STJ, AAO = r$AAO,
       bmiz = bmiz1(r$ht, r$wt, r$age_y, r$sex), z = as.list(z), d0 = as.list(d0), d2 = as.list(d2))
})
writeLines(paste0("window.LS_TESTS = ", toJSON(tests, auto_unbox = TRUE, digits = NA), ";"),
           file.path(out_dir, "ls_reference_cases.js"))
cat(sprintf("%d reference cases written (z, median and z = +2 diameter from z_score()).\n", length(tests)))
cat("Copy lauzon_schnittka.json/.js to models/ and ls_reference_cases.js to tests/.\n")
