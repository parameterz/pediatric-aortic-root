/* Calculator page. Needs models/root_bccg.js, tests/root_bccg_testcases.js, js/aortic-model.js. */
(function () {
  "use strict";
  const { M, R, KEYS, INFO, PHN, MAH, CSA, BANDS, interp, wfhOf, ref, zOf, dAt, medMm, haycock,
          category, fmtZ, mm, pctOf, sgnPct, el, selfTest } = AorticModel;

  // ---------- inputs ----------
  const $ = (id) => document.getElementById(id);
  const num = (id) => { const v = parseFloat(String($(id).value).replace(",", ".")); return isFinite(v) ? v : NaN; };
  function readAge() {
    if ($("am-d").checked) {
      const dob = $("dob").value, exam = $("exam").value;
      if (!dob || !exam) return NaN;
      return (new Date(exam + "T00:00:00Z") - new Date(dob + "T00:00:00Z")) / 86400000 / 365.25;
    }
    return num("age");
  }
  function read() {
    const d = {}; KEYS.forEach((k) => (d[k] = num("d-" + k)));
    return { d, ht: num("ht"), wt: num("wt"), age: readAge(), male: $("sex-m").checked };
  }

  // ---------- drawing ----------
  function drawScale(svg, r, z, dMm, color) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = 440, x0 = 14, x1 = W - 14, zMin = -3, zMax = 6, top = 26, h = 18;
    const X = (zz) => x0 + (Math.min(Math.max(zz, zMin), zMax) - zMin) / (zMax - zMin) * (x1 - x0);
    BANDS.forEach(([a, b, c]) => svg.appendChild(el("rect", { x: X(a), y: top, width: X(b) - X(a), height: h, style: `fill:var(${c})` })));
    svg.appendChild(el("rect", { x: x0, y: top, width: x1 - x0, height: h, style: "fill:none;stroke:var(--line)" }));
    [-2, 0, 2, 3.5, 5].forEach((t) => {
      const x = X(t);
      svg.appendChild(el("line", { x1: x, x2: x, y1: top, y2: top + h + 4, style: "stroke:var(--muted);stroke-width:1" }));
      svg.appendChild(el("text", { x, y: top + h + 16, "text-anchor": "middle", style: "fill:var(--ink);font:500 10.5px var(--font-num)" },
        r ? dAt(r, t).toFixed(1) : ""));
      svg.appendChild(el("text", { x, y: top + h + 28, "text-anchor": "middle", style: "fill:var(--muted);font:9.5px var(--font-num)" },
        "z " + (t === 0 ? "0" : (t > 0 ? "+" : "−") + Math.abs(t))));
    });
    if (isFinite(z)) {
      const x = X(z), off = z > zMax || z < zMin;
      svg.appendChild(el("line", { x1: x, x2: x, y1: top - 5, y2: top + h, style: `stroke:var(${color});stroke-width:2.5` }));
      svg.appendChild(el("circle", { cx: x, cy: top - 8, r: 3.5, style: `fill:var(${color})` }));
      const anchor = x > W - 90 ? "end" : x < 90 ? "start" : "middle";
      svg.appendChild(el("text", { x, y: top - 14, "text-anchor": anchor, style: `fill:var(${color});font:600 11px var(--font-num)` },
        dMm + " mm" + (off ? " (off scale)" : "")));
    }
  }

  // ---------- results ----------
  const rows = {};
  KEYS.forEach((k) => {
    const row = document.createElement("div"); row.className = "rrow";
    row.innerHTML =
      `<div class="rhead"><div class="rname"><span class="dot" style="background:var(${INFO[k].color})"></span>${INFO[k].name}<small>${INFO[k].short}</small></div>` +
      `<div class="zgroup"><span class="zval"><small>z</small><span class="zt">—</span></span><span class="chip none">Not entered</span></div></div>` +
      `<div class="rfacts"><span>Expected <b class="exp">—</b></span><span>Normal range <b class="rng">—</b></span><span>Percent of expected <b class="ratio">—</b></span></div>` +
      `<div class="rnote" hidden></div>` +
      `<svg viewBox="0 0 440 76" role="img" aria-label="${INFO[k].name} position on the normal range"></svg>`;
    $("results").appendChild(row);
    rows[k] = { row, zt: row.querySelector(".zt"), chip: row.querySelector(".chip"), exp: row.querySelector(".exp"), rng: row.querySelector(".rng"), ratio: row.querySelector(".ratio"), note: row.querySelector(".rnote"), svg: row.querySelector("svg") };
  });

  // ---------- comparison ----------
  function compareTable(key, p, hasWt) {
    const dMm = p.d[key], dcm = dMm / 10, ht = p.ht / 100, sexKey = p.male ? "M" : "F";
    const out = [];
    const add = (name, src, z, cut, extra) => out.push(Object.assign({ name, src, z, cut }, extra || {}));
    const na = (name, src, note) => out.push({ name, src, z: NaN, cut: NaN, note });

    if (hasWt) { const r = ref(key, p, true); add("This model", "Height, weight for height, age, sex; fitted shape", zOf(dMm, r), dAt(r, 2), { primary: true }); }
    else na("This model", "Height, weight for height, age, sex; fitted shape", "needs weight");
    const rc = ref(key, p, false);
    add("Height and sex", "Same model without weight", zOf(dMm, rc), dAt(rc, 2), { primary: !hasWt });
    if (PHN[key]) {
      const src = "PHN Echo Z-Score Project table (Lopez 2017): diameter ÷ BSA^0.5, Haycock BSA";
      if (!hasWt) na("PHN z-score", src, "needs weight");
      else { const rb = Math.sqrt(haycock(p.ht, p.wt)), [m, s] = PHN[key]; add("PHN z-score", src, (dcm / rb - m) / s, (m + 2 * s) * rb * 10); }
    }

    if (MAH[key]) {
      const src = "PHN values from Mahgerefteh 2021";
      if (!hasWt) na("BSA index, ≥ 2 y", "Diameter ÷ BSA^0.5 (Haycock); " + src, "needs weight");
      else if (p.age < 2) na("BSA index, ≥ 2 y", "Diameter ÷ BSA^0.5 (Haycock); " + src, "published for ≥ 2 y");
      else { const rb = Math.sqrt(haycock(p.ht, p.wt)), [m, s] = MAH[key].bsa; add("BSA index, ≥ 2 y", "Diameter ÷ BSA^0.5 (Haycock); " + src, (dcm / rb - m) / s, (m + 2 * s) * rb * 10); }
      if (p.age < 2) na("Height index, ≥ 2 y", "Diameter ÷ height; " + src, "published for ≥ 2 y");
      else { const [m, s] = MAH[key].ht; add("Height index, ≥ 2 y", "Diameter ÷ height; " + src, (dcm / ht - m) / s, (m + 2 * s) * ht * 10); }
    }
    if (CSA[key]) {
      const c = CSA[key];
      const bin = c.bins.find((b) => p.age >= b.lo && p.age < b.hi) || c.bins[c.bins.length - 1];
      const [bm, bsd] = bin[sexKey];
      const csa = Math.PI * (dcm / 2) ** 2 / ht;
      const toMm = (v) => 2 * Math.sqrt(v * ht / Math.PI) * 10;
      const reg = (k) => k.int + k.age * p.age + (p.male ? k.male + k.ma * p.age : 0);
      add("CSA/ht, age-group table", "Chaudhary et al., Pediatr Cardiol 2026", (csa - bm) / bsd, toMm(bm + 2 * bsd));
      const ms = reg(c.sheet), mr = reg(c.refit);
      add("CSA/ht, continuous age", "Coefficients as distributed in the 2026 supplement", (csa - ms) / bsd, toMm(ms + 2 * bsd));
      add("CSA/ht, continuous age", "Coefficients refit to the public PHN data", (csa - mr) / bsd, toMm(mr + 2 * bsd));
    }
    return out;
  }

  function renderCompare(p, hasWt) {
    const host = $("cmp-body"); host.textContent = "";
    const ready = p.ht > 0 && p.age >= 0;
    let any = false;
    KEYS.forEach((key) => {
      if (!(ready && p.d[key] > 0)) return;
      any = true;
      const h = document.createElement("h3");
      h.innerHTML = `<span class="dot" style="background:var(${INFO[key].color})"></span>${INFO[key].name}, ${p.d[key]} mm`;
      host.appendChild(h);
      const wrap = document.createElement("div"); wrap.className = "table-wrap";
      const tbl = document.createElement("table");
      tbl.innerHTML = `<thead><tr><th>Method</th><th class="num">z</th><th class="num">Diameter at z = +2</th></tr></thead>`;
      const tb = document.createElement("tbody");
      compareTable(key, p, hasWt).forEach((r) => {
        const tr = document.createElement("tr"); if (r.primary) tr.className = "primary";
        const td1 = document.createElement("td"); td1.textContent = r.name;
        const s = document.createElement("span"); s.className = "src"; s.textContent = r.src; td1.appendChild(s);
        const td2 = document.createElement("td"); td2.className = "num";
        const td3 = document.createElement("td"); td3.className = "num";
        if (isFinite(r.z)) { td2.textContent = fmtZ(r.z); td3.textContent = mm(r.cut); }
        else { td2.textContent = "—"; td3.textContent = r.note; td3.classList.add("sub"); }
        tr.append(td1, td2, td3); tb.appendChild(tr);
      });
      tbl.appendChild(tb); wrap.appendChild(tbl); host.appendChild(wrap);
    });
    if (!any) { const pEl = document.createElement("p"); pEl.className = "sub"; pEl.textContent = "Enter at least one diameter with height and age to compare."; host.appendChild(pEl); }
  }

  // ---------- weight for height ----------
  $("wf-rng").textContent = sgnPct(pctOf(R.wfh[0]), 0) + " to " + sgnPct(pctOf(R.wfh[1]), 0);
  $("wfh-legend").innerHTML = KEYS.map((k) => `<span><span class="dot" style="background:var(${INFO[k].color})"></span>${INFO[k].short}</span>`).join("");

  function drawWfh(childW) {
    const svg = $("wfh-svg");
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = 440, H = 220, l = 46, r = 14, t = 14, b = 40;
    const xMin = -50, xMax = 50, yMin = -12, yMax = 4;
    const X = (v) => l + (v - xMin) / (xMax - xMin) * (W - l - r);
    const Y = (v) => t + (yMax - v) / (yMax - yMin) * (H - t - b);
    for (let v = yMin; v <= yMax; v += 4) {
      svg.appendChild(el("line", { x1: l, x2: W - r, y1: Y(v), y2: Y(v), style: `stroke:var(--line);stroke-width:${v === 0 ? 1.5 : 1}` }));
      svg.appendChild(el("text", { x: l - 6, y: Y(v) + 4, "text-anchor": "end", style: "fill:var(--muted);font:10px var(--font-num)" }, (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) + "%"));
    }
    for (let v = xMin; v <= xMax; v += 25) {
      svg.appendChild(el("line", { x1: X(v), x2: X(v), y1: H - b, y2: H - b + 5, style: "stroke:var(--muted)" }));
      svg.appendChild(el("text", { x: X(v), y: H - b + 18, "text-anchor": "middle", style: "fill:var(--muted);font:10px var(--font-num)" }, (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) + "%"));
    }
    svg.appendChild(el("text", { x: (l + W - r) / 2, y: H - 4, "text-anchor": "middle", style: "fill:var(--muted);font:11px var(--font-ui)" }, "Weight for height"));
    KEYS.forEach((k) => {
      const f = M.structures[k].full.mu_wfh; let d = "";
      M.wfh.forEach((w, i) => { const px = pctOf(w); if (px < xMin || px > xMax) return; d += (d ? "L" : "M") + X(px).toFixed(1) + "," + Y(pctOf(f[i])).toFixed(1); });
      svg.appendChild(el("path", { d, style: `fill:none;stroke:var(${INFO[k].color});stroke-width:2` }));
    });
    if (childW != null && isFinite(childW)) {
      const cw = Math.min(Math.max(childW, R.wfh[0]), R.wfh[1]);
      const px = Math.min(Math.max(pctOf(cw), xMin), xMax);
      svg.appendChild(el("line", { x1: X(px), x2: X(px), y1: t, y2: H - b, style: "stroke:var(--ink);stroke-dasharray:3 3" }));
      KEYS.forEach((k) => svg.appendChild(el("circle", { cx: X(px), cy: Y(pctOf(interp(cw, M.wfh, M.structures[k].full.mu_wfh))), r: 4, style: `fill:var(${INFO[k].color})` })));
      const anchor = X(px) > W - 90 ? "end" : "start";
      svg.appendChild(el("text", { x: X(px) + (anchor === "start" ? 8 : -8), y: t + 10, "text-anchor": anchor, style: "fill:var(--ink);font:600 11px var(--font-num)" }, "this child"));
    }
  }

  function renderWfh(p, hasWt) {
    const tb = $("wf-tab"); tb.textContent = "";
    if (!(p.ht > 0)) { ["wf-typ", "wf-ent", "wf-pct"].forEach((id) => ($(id).textContent = "—")); drawWfh(null); return; }
    const L = Math.log(p.ht / 100);
    $("wf-typ").textContent = Math.exp(interp(L, M.lht, M.elw)).toFixed(1) + " kg";
    if (!hasWt) {
      $("wf-ent").textContent = "not entered"; $("wf-pct").textContent = "—";
      const tr = document.createElement("tr"); const td = document.createElement("td"); td.colSpan = 3; td.className = "sub";
      td.textContent = "Enter a weight to see its effect."; tr.appendChild(td); tb.appendChild(tr);
      drawWfh(null); return;
    }
    const w = wfhOf(p.ht, p.wt);
    $("wf-ent").textContent = p.wt.toFixed(1) + " kg";
    $("wf-pct").textContent = sgnPct(pctOf(w), 0) + " vs typical";
    KEYS.forEach((k) => {
      const f = M.structures[k].full, adj = interp(w, M.wfh, f.mu_wfh), sd = interp(L, M.lht, f.sigma);
      const tr = document.createElement("tr");
      tr.innerHTML = `<td><span class="dot" style="background:var(${INFO[k].color})"></span>${INFO[k].name}</td><td class="num">${sgnPct(pctOf(adj), 1)}</td><td class="num">${fmtZ(-adj / sd)}</td>`;
      tb.appendChild(tr);
    });
    drawWfh(w);
  }

  // ---------- main update ----------
  function update() {
    const p = read(), flags = [];
    $("age-out").textContent = $("am-d").checked && isFinite(p.age) ? "Age at exam: " + p.age.toFixed(2) + " years" : "";
    const ready = p.ht > 0 && isFinite(p.age) && p.age >= 0;
    const hasWt = p.wt > 0;

    if (ready) {
      if (p.ht < R.height_cm[0] || p.ht > R.height_cm[1]) flags.push(`Height is outside the reference range (${R.height_cm[0]} to ${R.height_cm[1]} cm), so results are extrapolations.`);
      if (p.age > R.age[1]) flags.push("Age is above 18 years, the upper limit of the reference group.");
      if (hasWt) { const w = wfhOf(p.ht, p.wt); if (w < R.wfh[0] || w > R.wfh[1]) flags.push("Weight for height is outside the reference range. Check the height and weight; results are extrapolations."); }
      KEYS.forEach((k) => { const v = p.d[k]; if (v > 0 && (v < 3 || v > 60)) flags.push(`The ${INFO[k].name.toLowerCase()} value looks unusual for millimetres. Check the units.`); });
      $("modelline").textContent = "Model: " + (hasWt ? "height, weight for height, age and sex." : "height and sex (enter a weight for the full model).");
    } else {
      $("modelline").textContent = "Enter height and age to see expected values.";
    }

    KEYS.forEach((k) => {
      const o = rows[k], dMm = p.d[k];
      if (!ready) {
        o.zt.textContent = "—"; o.chip.className = "chip none"; o.chip.textContent = "Not entered";
        o.exp.textContent = "—"; o.rng.textContent = "—"; o.ratio.textContent = "—"; o.note.hidden = true; o.row.classList.add("empty");
        drawScale(o.svg, null, NaN, NaN, INFO[k].color); return;
      }
      const r = ref(k, p, hasWt);
      const z = dMm > 0 ? zOf(dMm, r) : NaN;
      const [cls, label] = category(z);
      o.zt.textContent = isFinite(z) ? fmtZ(z) : "—";
      o.chip.className = "chip " + cls; o.chip.textContent = label;
      o.row.classList.toggle("empty", !isFinite(z));
      o.exp.textContent = mm(medMm(r));
      o.rng.textContent = dAt(r, -2).toFixed(1) + " to " + mm(dAt(r, 2));
      if (isFinite(z)) { const q = dMm / medMm(r); o.ratio.textContent = Math.round(100 * q) + "%"; } else o.ratio.textContent = "—";
      if (isFinite(z) && Math.abs(z) > 4) {
        o.note.hidden = false;
        o.note.textContent = `Beyond the reference data. Almost no healthy child in the reference group falls past z ±4, so this z-score extends the fitted shape rather than the data. The percent of expected (${Math.round(100 * dMm / medMm(r))}%) is the steadier number to follow.`;
      } else o.note.hidden = true;
      drawScale(o.svg, r, z, dMm, INFO[k].color);
    });

    const ul = $("flags"); ul.textContent = "";
    flags.forEach((f) => { const li = document.createElement("li"); li.textContent = f; ul.appendChild(li); });
    renderCompare(p, hasWt);
    renderWfh(p, hasWt);
  }

  // ---------- wiring ----------
  $("f").addEventListener("input", () => { $("example-note").hidden = true; update(); });
  $("f").addEventListener("submit", (e) => e.preventDefault());
  document.querySelectorAll('input[name="agemode"]').forEach((r) => r.addEventListener("change", () => {
    const dates = $("am-d").checked; $("age-years").hidden = dates; $("age-dates").hidden = !dates; update();
  }));
  $("clear").addEventListener("click", () => {
    ["d-SOV", "d-STJ", "d-AAO", "ht", "wt", "age", "dob", "exam"].forEach((id) => ($(id).value = ""));
    $("example-note").hidden = true; update(); $("d-SOV").focus();
  });

  // ---------- notes + self-test ----------
  $("formula").textContent = "z = ((D ÷ M)^λ − 1) ÷ (λ × σ)\nD at a given z = M × (1 + λσz)^(1/λ)";
  const mt = $("model-tab");
  KEYS.forEach((k) => {
    const st = M.structures[k], sg = st.full.sigma, tr = document.createElement("tr");
    tr.innerHTML = `<td><span class="dot" style="background:var(${INFO[k].color})"></span>${INFO[k].name}</td>` +
      `<td class="num">${st.n.toLocaleString("en-US")}</td><td class="num">${st.full.nu.toFixed(2)}</td>` +
      `<td class="num">${(100 * Math.min(...sg)).toFixed(1)}–${(100 * Math.max(...sg)).toFixed(1)}%</td>` +
      `<td class="num">${sgnPct(100 * (Math.exp(st.full.male) - 1), 1)}</td>`;
    mt.appendChild(tr);
  });
  const st = $("selftest"), res = selfTest(); st.textContent = "";
  const b = document.createElement("span"); b.className = res.pass ? "ok" : "bad";
  b.textContent = res.pass ? `${2 * res.n} of ${2 * res.n} reference z-scores match the R models` : "Mismatch with the R models";
  const a = document.createElement("a"); a.href = "validation.html"; a.textContent = "See every case";
  st.append(b, document.createTextNode(` (${res.n} children, full and height-and-sex versions; largest difference in z: ${res.maxd.toExponential(1)}). `), a, document.createTextNode("."));

  update();
})();
