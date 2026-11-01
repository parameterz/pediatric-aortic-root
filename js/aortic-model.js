/* Shared engine for the Pediatric Aortic Z-Scores pages.
 * Needs models/root_bccg.js (window.AORTIC_MODEL) loaded first; tests/root_bccg_testcases.js
 * (window.AORTIC_TESTS) is optional. Classic script (no modules) so pages also open from disk.
 *
 * Model per level (BCCG / LMS-type GAMLSS):
 *   ln M = f(ln height) + g(weight for height) + h(age) + sex     (full model)
 *   ln M = f(ln height) + sex                                     (height-and-sex model)
 *   sigma = s(ln height);  nu (lambda) constant
 *   z = ((D/M)^nu - 1) / (nu * sigma);  D at z = M * (1 + nu*sigma*z)^(1/nu)
 */
(function (global) {
  "use strict";
  const M = global.AORTIC_MODEL;
  if (!M) throw new Error("AORTIC_MODEL not loaded: include models/root_bccg.js before js/aortic-model.js");
  const TESTS = global.AORTIC_TESTS || [];
  const R = M.meta.ranges;
  const KEYS = ["SOV", "STJ", "AAO"];
  const INFO = {
    SOV: { name: "Sinuses of Valsalva", short: "SOV", color: "--c-sov" },
    STJ: { name: "Sinotubular junction", short: "STJ", color: "--c-stj" },
    AAO: { name: "Ascending aorta", short: "AAo", color: "--c-aao" }
  };

  // ---- comparison methods (published values) ----
  // PHN Echo Z-Score Project table: D / BSA^0.5 (Haycock), [mean, SD], all ages
  const PHN = { SOV: [2.06, 0.18], STJ: [1.69, 0.16], AAO: [1.79, 0.18] };
  // Mahgerefteh 2021 Table 2 (>= 2 y): [mean, SD] of D/Ht and D/BSA^0.5
  const MAH = { SOV: { ht: [1.61, 0.15], bsa: [2.07, 0.18] }, AAO: { ht: [1.39, 0.14], bsa: [1.78, 0.18] } };
  // Chaudhary 2026 CSA/ht age-group tables [mean, SD] and regression coefficients
  const CSA = {
    SOV: {
      bins: [
        { lo: 0, hi: 1 / 12, F: [1.4128447, 0.2550527], M: [1.5346930, 0.2775356] },
        { lo: 1 / 12, hi: 3, F: [1.8518291, 0.4176144], M: [1.9816039, 0.4120444] },
        { lo: 3, hi: 6, F: [2.2326539, 0.3698699], M: [2.4048436, 0.3998567] },
        { lo: 6, hi: 12, F: [2.5486409, 0.4543347], M: [2.8695985, 0.5095964] },
        { lo: 12, hi: 16, F: [2.9936960, 0.4995857], M: [3.2347010, 0.5950158] },
        { lo: 16, hi: Infinity, F: [3.1964701, 0.5267062], M: [3.6397332, 0.6259371] }
      ],
      sheet: { int: 1.6414159903, age: 0.1294167470, male: 0.0961634070, ma: 0.0147000314 },
      refit: { int: 1.6413, age: 0.0962, male: 0.1302, ma: 0.0146 }
    },
    AAO: {
      bins: [
        { lo: 0, hi: 1 / 12, F: [1.1566680, 0.2718176], M: [1.2500483, 0.2862416] },
        { lo: 1 / 12, hi: 3, F: [1.3894141, 0.3164411], M: [1.3956703, 0.3201912] },
        { lo: 3, hi: 6, F: [1.6674859, 0.3513082], M: [1.7176683, 0.3666948] },
        { lo: 6, hi: 12, F: [1.9498509, 0.3534838], M: [2.0984815, 0.4530524] },
        { lo: 12, hi: 16, F: [2.3370879, 0.4802758], M: [2.3480720, 0.4751761] },
        { lo: 16, hi: Infinity, F: [2.4963250, 0.5387956], M: [2.6417850, 0.5027709] }
      ],
      sheet: { int: 1.2612332551, age: 0.0538404716, male: 0.0751494789, ma: 0.0031328790 },
      refit: { int: 1.2612, age: 0.0751, male: 0.0553, ma: 0.0031 }
    }
  };

  // ---- math ----
  function interp(x, xs, ys) {
    const n = xs.length;
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
    const t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
  }
  const wfhOf = (htCm, wtKg) => Math.log(wtKg) - interp(Math.log(htCm / 100), M.lht, M.elw);
  const typicalWeight = (htCm) => Math.exp(interp(Math.log(htCm / 100), M.lht, M.elw));

  // Reference for one child: lm = ln(median, cm), sigma, nu. p = {ht (cm), wt (kg), age (y), male (bool)}
  function ref(key, p, hasWt) {
    const s = M.structures[key], L = Math.log(p.ht / 100);
    if (hasWt) {
      const f = s.full, wfh = wfhOf(p.ht, p.wt);
      const lm = interp(L, M.lht, f.mu_lht) + interp(wfh, M.wfh, f.mu_wfh) +
                 interp(p.age, M.age, f.mu_age) + (p.male ? f.male : 0);
      return { lm, sigma: interp(L, M.lht, f.sigma), nu: f.nu, wfh };
    }
    const h = s.hs;
    return { lm: interp(L, M.lht, h.mu_lht) + (p.male ? h.male : 0), sigma: interp(L, M.lht, h.sigma), nu: h.nu, wfh: null };
  }
  // z for a diameter in mm under a reference r; optional nu override (shape explainer)
  function zOf(dMm, r, nu) {
    const v = nu == null ? r.nu : nu;
    const ratio = dMm / 10 / Math.exp(r.lm);
    return Math.abs(v) < 1e-6 ? Math.log(ratio) / r.sigma : (Math.pow(ratio, v) - 1) / (v * r.sigma);
  }
  // Diameter (mm) at a given z
  function dAt(r, z, nu) {
    const v = nu == null ? r.nu : nu;
    const med = 10 * Math.exp(r.lm);
    if (Math.abs(v) < 1e-6) return med * Math.exp(r.sigma * z);
    const base = 1 + v * r.sigma * z;
    return base > 0 ? med * Math.pow(base, 1 / v) : 0;
  }
  const medMm = (r) => 10 * Math.exp(r.lm);
  const haycock = (htCm, wtKg) => 0.024265 * Math.pow(htCm, 0.3964) * Math.pow(wtKg, 0.5378);

  // ---- presentation ----
  const BANDS = [[-3, -2, "--band-small"], [-2, 2, "--band-normal"], [2, 3.5, "--band-mild"], [3.5, 5, "--band-moderate"], [5, 6, "--band-severe"]];
  function category(z) {
    if (!isFinite(z)) return ["none", "Not entered"];
    if (z < -2) return ["small", "Below the normal range"];
    if (z < 2) return ["normal", "Within the normal range"];
    if (z < 3.5) return ["mild", "Mild dilation"];
    if (z < 5) return ["moderate", "Moderate dilation"];
    return ["severe", "Severe dilation"];
  }
  const fmtZ = (z) => (z >= 0 ? "+" : "−") + Math.abs(z).toFixed(2);
  const mm = (v) => v.toFixed(1) + " mm";
  const pctOf = (w) => 100 * (Math.exp(w) - 1);
  const sgnPct = (v, d) => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(d) + "%";

  const SVGNS = "http://www.w3.org/2000/svg";
  function el(tag, attrs, text) {
    const e = document.createElementNS(SVGNS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (text != null) e.textContent = text;
    return e;
  }

  // Run the reference cases; returns { n, maxd, rows }
  function selfTest() {
    let maxd = 0;
    const rows = TESTS.map((t) => {
      const p = { ht: t.height_cm, wt: t.weight_kg, age: t.age_y, male: t.sex === "M" };
      const zf = zOf(t.diam_mm, ref(t.structure, p, true));
      const zh = zOf(t.diam_mm, ref(t.structure, p, false));
      const df = zf - t.z_R, dh = zh - t.z_R_hs;
      maxd = Math.max(maxd, Math.abs(df), Math.abs(dh));
      return Object.assign({}, t, { z_js: zf, z_js_hs: zh, diff: df, diff_hs: dh });
    });
    return { n: TESTS.length, maxd, rows, pass: maxd < 0.001 };
  }

  // ---- units for height and weight (stored values are always cm and kg) ----
  const CONV = { ht: { cm: 1, in: 2.54 }, wt: { kg: 1, lb: 0.45359237 } };
  const METRIC = { ht: "cm", wt: "kg" }, OTHER = { ht: "in", wt: "lb" };
  const UNITS_KEY = "aortic-units-v1";
  const units = { ht: "cm", wt: "kg" };
  try {
    const u = JSON.parse(localStorage.getItem(UNITS_KEY));
    if (u && CONV.ht[u.ht]) units.ht = u.ht;
    if (u && CONV.wt[u.wt]) units.wt = u.wt;
  } catch (e) { /* storage unavailable: metric */ }
  const parseNum = (v) => { const x = parseFloat(String(v == null ? "" : v).replace(",", ".")); return isFinite(x) ? x : NaN; };
  const toMetric = (kind, v) => v * CONV[kind][units[kind]];
  const fromMetric = (kind, v) => v / CONV[kind][units[kind]];
  function setUnit(kind, u) {
    units[kind] = u;
    try { localStorage.setItem(UNITS_KEY, JSON.stringify(units)); } catch (e) { /* ignore */ }
  }
  // Wire a field: <input id="X"> plus <button id="X-unit"> and optional <span id="X-alt">.
  // The input holds a metric example value in the HTML; it is converted if the saved preference is US units.
  function bindUnitField(inputId, kind, onChange) {
    const inp = document.getElementById(inputId), btn = document.getElementById(inputId + "-unit"), alt = document.getElementById(inputId + "-alt");
    if (!inp || !btn) return;
    const name = kind === "ht" ? "Height" : "Weight";
    if (units[kind] !== METRIC[kind]) { const v = parseNum(inp.value); if (isFinite(v)) inp.value = (v / CONV[kind][units[kind]]).toFixed(1); }
    const paint = () => {
      const next = units[kind] === METRIC[kind] ? OTHER[kind] : METRIC[kind];
      btn.textContent = units[kind];
      btn.setAttribute("aria-label", `${name} in ${units[kind]}. Switch to ${next}`);
      btn.title = `Switch to ${next}`;
      if (alt) { const v = parseNum(inp.value); alt.textContent = units[kind] !== METRIC[kind] && isFinite(v) ? `= ${toMetric(kind, v).toFixed(1)} ${METRIC[kind]}` : ""; }
    };
    btn.addEventListener("click", () => {
      const v = parseNum(inp.value), m = isFinite(v) ? toMetric(kind, v) : NaN;
      setUnit(kind, units[kind] === METRIC[kind] ? OTHER[kind] : METRIC[kind]);
      if (isFinite(m)) inp.value = fromMetric(kind, m).toFixed(1);
      paint(); if (onChange) onChange();
    });
    inp.addEventListener("input", paint);
    paint();
  }

  // ---- height-for-age check (needs models/age_height.js; silent if it is missing) ----
  const AH = global.AORTIC_AGEHT || null;
  function htForAgeZ(htCm, age, male) {
    if (!AH) return NaN;
    const s = AH[male ? "M" : "F"];
    const mu = interp(age, AH.age, s.mu), sg = interp(age, AH.age, s.sigma), nu = s.nu;
    return Math.abs(nu) < 1e-6 ? Math.log(htCm / mu) / sg : (Math.pow(htCm / mu, nu) - 1) / (nu * sg);
  }
  function ageCheck(p) {
    if (!AH || !(p.ht > 0) || !(p.age >= 0) || p.age > 18) return null;
    const z = htForAgeZ(p.ht, p.age, p.male);
    if (!(Math.abs(z) > 4)) return null;
    const typical = interp(p.age, AH.age, AH[p.male ? "M" : "F"].mu);
    let msg = `A height of ${p.ht.toFixed(0)} cm is unusual for age ${+p.age.toFixed(1)} (typical about ${typical.toFixed(0)} cm). Check the age and height.`;
    if (p.age >= 1 && Math.abs(htForAgeZ(p.ht, p.age / 12, p.male)) < 2.5) msg += ` If the age was entered in months, it is ${(p.age / 12).toFixed(1)} years.`;
    return msg + " Very tall or short stature, as in Marfan or Turner syndrome, can also cause this.";
  }

  // ---- site version and author (one place to update; added to every page's footer) ----
  const SITE = {
    version: "0.2", status: "early version", date: "31 October 2026",
    author: "Dan Dyar, MA, ACS, RDCS, FASE",
    link: "https://www.linkedin.com/in/dan-dyar-ma-acs-rdcs-ae-pe-fe-fase/"
  };
  function addFooter() {
    const f = document.querySelector("footer");
    if (!f || f.querySelector(".site-credit")) return;
    const p = document.createElement("p");
    p.className = "site-credit";
    const a = document.createElement("a");
    a.href = SITE.link; a.target = "_blank"; a.rel = "noopener"; a.textContent = SITE.author;
    p.append("Created by ", a, ` · Version ${SITE.version} (${SITE.status}), ${SITE.date}`);
    f.prepend(p);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", addFooter); else addFooter();

  global.AorticModel = {
    SITE, units, toMetric, fromMetric, setUnit, bindUnitField, parseNum, htForAgeZ, ageCheck,
    M, TESTS, R, KEYS, INFO, PHN, MAH, CSA, BANDS,
    interp, wfhOf, typicalWeight, ref, zOf, dAt, medMm, haycock,
    category, fmtZ, mm, pctOf, sgnPct, el, selfTest
  };
})(window);
