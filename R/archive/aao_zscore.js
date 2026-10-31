// AAo z-score from the PHN location-scale GAM, evaluated from lookup tables.
// Usage:
//   const model = await loadAaoModel("aao_model.json");
//   const r = aaoZ(model, { aadMm: 24, heightCm: 128, weightKg: 28, ageYears: 8, sex: "F" });

export async function loadAaoModel(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Could not load model: ${res.status}`);
  return res.json();
}

// Linear interpolation on an increasing grid; clamps outside the range (like R's rule = 2)
function interp(x, xs, ys) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid; else hi = mid;
  }
  const t = (x - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] + t * (ys[hi] - ys[lo]);
}

export function aaoZ(model, { aadMm, heightCm, weightKg, ageYears, sex }) {
  const flags = [];
  const r = model.meta.ranges;
  if (heightCm < r.height_cm[0] || heightCm > r.height_cm[1])
    flags.push(`Height outside reference range (${r.height_cm[0]}–${r.height_cm[1]} cm)`);
  if (ageYears < r.age[0] || ageYears > r.age[1])
    flags.push(`Age outside reference range (${r.age[0]}–${r.age[1]} y)`);

  const L = Math.log(heightCm / 100);
  const wfh = Math.log(weightKg) - interp(L, model.lht, model.elw);
  if (wfh < r.wfh[0] || wfh > r.wfh[1])
    flags.push("Weight-for-height outside reference range: check height and weight");

  const mu =
    interp(L, model.lht, model.f_lht) +
    interp(wfh, model.wfh, model.f_wfh) +
    interp(ageYears, model.age, model.f_age) +
    (sex === "M" ? model.male : 0);
  const sd = interp(L, model.lht, model.sd);
  const z = (Math.log(aadMm / 10) - mu) / sd;

  return {
    z,
    predictedMm: 10 * Math.exp(mu),
    lower2Mm: 10 * Math.exp(mu - 2 * sd),
    upper2Mm: 10 * Math.exp(mu + 2 * sd),
    weightForHeightPct: 100 * (Math.exp(wfh) - 1),   // % above/below typical weight for height
    flags
  };
}
