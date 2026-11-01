# Pediatric Aortic Z-Scores

*An analysis of the Pediatric Heart Network public database*

A browser-based z-score calculator for the **sinuses of Valsalva**, **sinotubular junction** and **ascending aorta** in children, built from the NIH/NHLBI Pediatric Heart Network (PHN) Echo-Z normal echocardiogram dataset (about 3,200 healthy, non-obese children from birth to 18 years).

**Site:** enable GitHub Pages for this repository, or open `index.html` directly from disk.

| Page | What it shows |
|---|---|
| `index.html` | Calculator: z-score, expected diameter, normal range and percent of expected for each level, with a weight-for-height panel and a comparison with other methods |
| `charts.html` | Reference charts of diameter by height (z −2 to +5) for each level and sex, plus a visit tracker |
| `shape.html` | Interactive explainer: how the Box-Cox shape λ changes the z-score of one diameter |
| `diagnostics.html` | Residuals against body size, Q-tests, Q-Q, histogram, worm plots, shape comparison and local height exponents (aggregated summaries only) |
| `validation.html` | Every reference case recomputed in the browser and compared with R |

Everything runs in the browser. Nothing entered on the pages is sent anywhere.

## What the calculator gives

For each level of the aorta:

- **z-score**, with dilation bands from the AHA scientific statement on aortopathy in children (normal < 2, mild 2 to < 3.5, moderate 3.5 to < 5, severe ≥ 5)
- **expected diameter** and the **normal range** (z −2 to +2) in mm
- **percent of expected** (100 × measured ÷ predicted median), which does not depend on the tail model and is the steadier number for severe dilation
- a note when |z| > 4, where almost no reference children exist and any z-score is an extrapolation
- a **weight-for-height** panel showing how build changes the expected diameter
- a collapsible **comparison** with other methods: the PHN BSA z-score (Lopez 2017), height and BSA indices (Mahgerefteh 2021), Lauzon-Schnittka (2025; height, weight, WHO BMI-for-age z and age, from a large clinical cohort that includes overweight and obese children), and CSA/height (Chaudhary 2026)

## The model

Each level has its own **LMS-type model**: GAMLSS with the Box-Cox Cole-Green (BCCG) distribution, fitted in R with the `gamlss` package.

| LMS part | Parameter | Model |
|---|---|---|
| **M**, the median | μ (log link) | smooth log height + smooth weight for height + smooth age + sex |
| **S**, the spread (≈ CV) | σ (log link) | smooth log height |
| **L**, the shape | ν (Box-Cox power, λ) | constant per level, estimated from the data |

z = ((D ÷ M)^λ − 1) ÷ (λ × σ), and the diameter at a given z is M × (1 + λσz)^(1/λ).

Weight for height is ln(weight) minus the smoothed average ln(weight) for that height in the same reference group. Without a weight, the calculator uses the same kind of model built from height and sex only.

| Level | Children | λ | σ | Male vs female |
|---|---|---|---|---|
| Sinuses of Valsalva | 3,210 | 0.05 (≈ log-normal) | 8.1–9.0% | +4.1% |
| Sinotubular junction | 3,209 | 0.39 | 8.3–10.2% | +2.1% |
| Ascending aorta | 3,207 | 0.37 | 9.4–10.9% | +1.5% |

Fixing λ at 2 (the scale of an area, as in CSA/height) was strongly rejected for all three levels.

### Checks

- z-scores have mean 0 and SD 1 in every height decile and age group, with no remaining trend against height, age, weight, BMI or sex
- 2.1–2.6% of children beyond ±2 on each side (2.3% expected)
- Royston–Wright Q-tests by height group: spread, skewness and kurtosis pass at every level, except slightly heavy tails at the sinotubular junction
- the browser implementation reproduces the R z-scores for 60 reference cases to within 1e-4 (self-test shown on the page)

### Data cleaning

Only clear errors were excluded, never values that were merely unusual.

1. **Body size:** five records with physiologically impossible height–weight combinations were excluded from all models.
2. **Measurements:** a reading was excluded only when it was implausible **and isolated to one level** while the child's other levels were normal: one sinotubular junction and two ascending aorta readings. A child large at all three levels was kept, since agreement across levels points to real anatomy.

## Repository layout

```
index.html, charts.html, shape.html, diagnostics.html, validation.html
css/site.css                     shared styles
js/aortic-model.js               shared engine: tables -> median, sigma, lambda -> z and cutoffs
js/<page>.js                     one script per page
models/root_bccg.json / .js      exported model tables (the .js wrapper is what the pages load)
models/diagnostics.json / .js    aggregated diagnostics (from R/export_diagnostics.R)
tests/root_bccg_testcases.csv / .js   reference cases for the self-test and validation page
R/aortic_root_models.R           data preparation, cleaning, exponents, location-scale GAMs, checks
R/export_bccg.R                  fits the BCCG models and exports the tables and test cases
R/export_diagnostics.R           aggregated diagnostics for diagnostics.html
R/export_age_height.R            height-for-age reference used to flag an implausible age (e.g., months entered as years)
models/age_height.json / .js     output of export_age_height.R
models/lauzon_schnittka.js       Lauzon-Schnittka 2025 coefficients, spline tables and WHO BMI-for-age LMS tables, from the authors' published Excel calculator
js/lauzon-schnittka.js           their z-score formulas; checked against the workbook (tests/ls_reference_cases.js)
R/diagnostics_qtests.R           console version of the shape comparison, Q-tests and worm plots
R/archive/                       earlier ascending-aorta exploration (CSA/height critique, first calculator)
```

The pages load data through classic `<script>` tags (`window.AORTIC_MODEL`, `window.AORTIC_TESTS`, `window.AORTIC_DIAG`) rather than `fetch`, so they also work when opened straight from disk.

## Reproducing the models

The PHN data are **not** included in this repository. Download the Echo-Z/Normal ECG public use dataset from [pediatricheartnetwork.org](https://www.pediatricheartnetwork.org), load it in R as `echo_zscores_phn`, then run:

```r
source("R/aortic_root_models.R")    # needs mgcv, jsonlite, ggplot2
source("R/export_bccg.R")           # needs gamlss
source("R/export_diagnostics.R")    # aggregated diagnostics for the web page
source("R/export_age_height.R")     # height-for-age curves for the age check
```

Outputs are written to `aao_explore/`. Copy `root_bccg.json` and `root_bccg.js` to `models/`, `root_bccg_testcases.csv` and `.js` to `tests/`, and `diagnostics.json`/`.js` and `age_height.json`/`.js` to `models/`.

Height and weight can be entered in cm/kg or in/lb (the small unit button inside each field); the choice is remembered in the browser, and all calculations use cm and kg.

## Limits

- The reference group excluded obesity; results far outside the reference weight-for-height range are extrapolations.
- Diameters were recorded to the nearest millimetre; in newborns, 1 mm moves z by about 1.
- Calibration is close to exact out to about z ±2.5. Beyond ±3 the tails hold roughly twice as many reference children as a normal distribution predicts, and past about ±4 any z-score extends the fitted shape rather than the data.
- For research and education. It does not replace clinical judgment.

## Acknowledgment

The NIH/NHLBI Pediatric Heart Network Echo-Z/Normal ECG dataset was used in preparation of this work. Data were downloaded from www.pediatricheartnetwork.org on 01/01/2026.

This is an independent analysis and is not a product of the Pediatric Heart Network.

## Key references

- Lopez L, et al. Circ Cardiovasc Imaging. 2017;10:e006979. doi:10.1161/CIRCIMAGING.117.006979
- Mahgerefteh J, et al. Pediatr Cardiol. 2021;42:1284–1292. doi:10.1007/s00246-021-02609-x
- Chaudhary AA, et al. Pediatr Cardiol. 2026. doi:10.1007/s00246-026-04525-4
- Cole TJ, Green PJ. Stat Med. 1992;11:1305–1319 (the LMS method)
- Rigby RA, Stasinopoulos DM. J R Stat Soc C. 2005;54:507–554 (GAMLSS)
- Lauzon-Schnittka J, et al. Circ Cardiovasc Imaging. 2025;18:e017944. doi:10.1161/CIRCIMAGING.124.017944
- Herrera-Leaño N, et al. J Am Heart Assoc. 2026;15:e050142 (SEAMLSS)
- Morris SA, et al. Circulation. 2024 (AHA scientific statement on aortopathy in children)

## License

MIT. See `LICENSE`.
