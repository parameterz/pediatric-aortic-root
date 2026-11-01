/* Lauzon-Schnittka et al. (Circ Cardiovasc Imaging 2025;18:e017944) z-scores.
 * A line-by-line port of the authors' R function z_score() (zscorecalculator.r).
 * Needs models/lauzon_schnittka.js (window.LS_MODEL), written by R/export_ls_model.R
 * straight from that function.
 *
 * Per level and sex:
 *   Ta     = age^p (ln age when p = 0), age floored at 0.001 y
 *   BMIz   = WHO BMI-for-age LMS z (daily table < 2 y, monthly table from 2 y), clamped to [-5, 5]
 *   mu     = exp(b1 + b2 ln(height cm) + b3 ln(weight kg) + b4 Ta + b5 BMIz + s1(Ta) + s2(BMIz))
 *   sigma  = exp(b6 + b7 Ta + s3(Ta));   lambda = b8 + b9 Ta + s4(Ta)
 *   z      = ((D / mu)^lambda - 1) / (lambda sigma)    (ln(D / mu) / sigma when lambda = 0)
 * s1..s4 are cubic B-splines (order 4) given by knots and weights; a missing spline is 0. */
(function (global) {
  "use strict";
  const LS = global.LS_MODEL;
  if (!LS) return;

  const lerp = (x, x1, y1, x2, y2) => y1 + (x - x1) * (y2 - y1) / (x2 - x1);

  // WHO BMI-for-age z, as bmiz_who() in the R code (unclamped). NaN outside the tables.
  function bmiz(male, age, htCm, wtKg) {
    const t = LS.bmiz[male ? "m" : "f"];
    let tab, i, x1, x2;
    if (age < 2) { tab = t.day; i = Math.floor(age * 365.25); x1 = i / 365.25; x2 = (i + 1) / 365.25; }
    else { tab = t.mon; const j = Math.floor(age * 12); i = j - 24; x1 = j / 12; x2 = (j + 1) / 12; }
    if (i < 0 || i + 1 >= tab.L.length) return NaN;
    const l = lerp(age, x1, tab.L[i], x2, tab.L[i + 1]);
    const m = lerp(age, x1, tab.M[i], x2, tab.M[i + 1]);
    const s = lerp(age, x1, tab.S[i], x2, tab.S[i + 1]);
    const bmi = wtKg / ((htCm / 100) * (htCm / 100));
    return l === 0 ? Math.log(bmi / m) / s : (Math.pow(bmi / m, l) - 1) / (l * s);
  }

  // B-spline basis, Cox-de Boor recursion exactly as Nik_x() in the R code (1-based i)
  function N(i, k, x, t) {
    if (k === 1) return x >= t[i - 1] && x <= t[i] ? 1 : 0;
    return N(i, k - 1, x, t) * (x - t[i - 1]) / (t[i + k - 2] - t[i - 1]) +
           N(i + 1, k - 1, x, t) * (t[i + k - 1] - x) / (t[i + k - 1] - t[i]);
  }
  function S(x, sp) {
    if (!sp || !sp.t || !sp.w) return 0;
    let s = 0;
    for (let i = 1; i <= sp.t.length - 4; i++) s += sp.w[i - 1] * N(i, 4, x, sp.t);
    return s;
  }

  // Reference for one child and level: { mu (mm), sigma, nu, bmiz }, or null outside the model's range
  function ref(key, p) {
    const m = LS.levels[key + (p.male ? "m" : "f")];
    if (!m || !(p.ht > 0) || !(p.wt > 0) || !(p.age >= 0)) return null;
    const bzRaw = bmiz(p.male, p.age, p.ht, p.wt);
    if (!isFinite(bzRaw)) return null;                       // beyond the WHO tables (about 19 y)
    const b = m.b, a = Math.max(p.age, 0.001);
    const Ta = m.p === 0 ? Math.log(a) : Math.pow(a, m.p);
    const bz = Math.min(5, Math.max(-5, bzRaw));
    const mu = Math.exp(b[0] + b[1] * Math.log(p.ht) + b[2] * Math.log(p.wt) + b[3] * Ta + b[4] * bz + S(Ta, m.s1) + S(bz, m.s2));
    const sigma = Math.exp(b[5] + b[6] * Ta + S(Ta, m.s3));
    const nu = b[7] + b[8] * Ta + S(Ta, m.s4);
    return { mu, sigma, nu, bmiz: bzRaw };
  }
  const z = (D, r) => r.nu === 0 ? Math.log(D / r.mu) / r.sigma : (Math.pow(D / r.mu, r.nu) - 1) / (r.nu * r.sigma);
  const dAt = (r, zz) => r.nu === 0 ? r.mu * Math.exp(r.sigma * zz) : r.mu * Math.pow(Math.max(0, 1 + r.nu * r.sigma * zz), 1 / r.nu);

  // Check against z_score() run in R (tests/ls_reference_cases.js): z, median and the z = +2 diameter
  function selfTest() {
    const T = global.LS_TESTS || [];
    let maxd = 0, maxb = 0, maxmm = 0;
    const rows = T.map((t) => {
      const p = { ht: t.ht, wt: t.wt, age: t.age_y, male: t.sex === "m" };
      const out = { t, z: {}, d2: {} };
      ["SOV", "STJ", "AAO"].forEach((k) => {
        const r = ref(k, p), zz = z(t[k], r);
        out.z[k] = zz; out.d2[k] = dAt(r, 2); out.bmiz = r.bmiz;
        maxd = Math.max(maxd, Math.abs(zz - t.z[k]));
        if (t.d0) maxmm = Math.max(maxmm, Math.abs(r.mu - t.d0[k]), Math.abs(out.d2[k] - t.d2[k]));
      });
      maxb = Math.max(maxb, Math.abs(out.bmiz - t.bmiz));
      return out;
    });
    return { n: T.length, rows, maxd, maxb, maxmm, pass: T.length > 0 && maxd < 1e-6 && maxb < 1e-6 && maxmm < 1e-6 };
  }

  global.LauzonSchnittka = { ref, z, dAt, bmiz, selfTest, meta: LS.meta };
})(window);
