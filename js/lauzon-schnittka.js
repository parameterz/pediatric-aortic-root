/* Lauzon-Schnittka et al. (Circ Cardiovasc Imaging 2025;18:e017944) z-scores, as implemented in the
 * authors' Excel calculator. Needs models/lauzon_schnittka.js (window.LS_MODEL).
 *
 * Per level and sex (sheet rows in the workbook):
 *   f(age) = age^p (ln age when p = 0), age floored at the first table age
 *   BMIz   = WHO BMI-for-age LMS z, clamped to [-5, 5] for the model
 *   mu     = exp(b1 + b2 ln(height cm) + b3 ln(weight kg) + b4 f(age) + b5 BMIz + s1(age) + s2(BMIz))
 *   sigma  = exp(b6 + b7 f(age) + s3(age));   nu = b8 + b9 f(age) + s4(age)
 *   z      = ((D / mu)^nu - 1) / (nu * sigma)
 * s1..s4 are spline tables, linearly interpolated as in the workbook. */
(function (global) {
  "use strict";
  const LS = global.LS_MODEL;
  if (!LS) return;
  function interp(x, xs, ys) {
    const n = xs.length;
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
    return ys[lo] + (x - xs[lo]) / (xs[hi] - xs[lo]) * (ys[hi] - ys[lo]);
  }
  // WHO BMI-for-age z (unclamped)
  function bmiz(male, age, htCm, wtKg) {
    const t = LS.bmiz[male ? "m" : "f"];
    const L = interp(age, t.age, t.L), Mv = interp(age, t.age, t.M), S = interp(age, t.age, t.S);
    const bmi = wtKg / Math.pow(htCm / 100, 2);
    return (Math.pow(bmi / Mv, L) - 1) / (L * S);
  }
  // Reference for one child and level: { mu (mm), sigma, nu, bmiz }
  function ref(key, p) {
    const m = LS.levels[key + (p.male ? "m" : "f")];
    if (!m || !(p.ht > 0) || !(p.wt > 0) || !(p.age >= 0)) return null;
    const b = m.b, a = Math.max(p.age, m.age[0]);
    const fa = m.p === 0 ? Math.log(a) : Math.pow(a, m.p);
    const bzRaw = bmiz(p.male, p.age, p.ht, p.wt), bz = Math.min(Math.max(bzRaw, -5), 5);
    const mu = Math.exp(b[0] + b[1] * Math.log(p.ht) + b[2] * Math.log(p.wt) + b[3] * fa + b[4] * bz +
                        interp(a, m.age, m.s1) + interp(bz, m.bmiz, m.s2));
    const sigma = Math.exp(b[5] + b[6] * fa + interp(a, m.age, m.s3));
    const nu = b[7] + b[8] * fa + interp(a, m.age, m.s4);
    return { mu, sigma, nu, bmiz: bzRaw };
  }
  const z = (D, r) => Math.abs(r.nu) < 1e-6 ? Math.log(D / r.mu) / r.sigma : (Math.pow(D / r.mu, r.nu) - 1) / (r.nu * r.sigma);
  const dAt = (r, zz) => Math.abs(r.nu) < 1e-6 ? r.mu * Math.exp(r.sigma * zz) : r.mu * Math.pow(Math.max(0, 1 + r.nu * r.sigma * zz), 1 / r.nu);

  function selfTest() {
    const T = global.LS_TESTS || [];
    let maxd = 0, maxb = 0;
    const rows = T.map((t) => {
      const p = { ht: t.ht, wt: t.wt, age: t.age_y, male: t.sex === "m" };
      const out = { t, z: {} };
      ["SOV", "STJ", "AAO"].forEach((k) => {
        const r = ref(k, p), zz = z(t[k], r);
        out.z[k] = zz; out.bmiz = r.bmiz;
        maxd = Math.max(maxd, Math.abs(zz - t.z[k]));
      });
      maxb = Math.max(maxb, Math.abs(out.bmiz - t.bmiz));
      return out;
    });
    return { n: T.length, rows, maxd, maxb, pass: maxd < 1e-6 && maxb < 1e-6 };
  }

  global.LauzonSchnittka = { ref, z, dAt, bmiz, selfTest, meta: LS.meta };
})(window);
