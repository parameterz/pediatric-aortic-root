/* Reference charts and visit tracker.
 * Curves: diameter vs height at fixed z, from the height-and-sex model (a child of typical build).
 * Visits: z and percent of expected from the full model when weight is entered. */
(function () {
  "use strict";
  const { R, KEYS, INFO, ref, zOf, dAt, medMm, category, fmtZ, el } = AorticModel;
  const $ = (id) => document.getElementById(id);
  const ZS = [-2, 0, 2, 3.5, 5];
  const FILLS = [[-3, -2, "--band-small"], [-2, 2, "--band-normal"], [2, 3.5, "--band-mild"], [3.5, 5, "--band-moderate"], [5, 6, "--band-severe"]];
  const H0 = Math.max(45, Math.ceil(R.height_cm[0])), H1 = Math.min(205, Math.floor(R.height_cm[1]));
  const state = { key: "AAO", male: true };

  // ---------- visits ----------
  const EXAMPLE = [
    { age: 6, ht: 116, wt: 21, SOV: 19, STJ: 15, AAO: 21 },
    { age: 8, ht: 128, wt: 28, SOV: 21, STJ: 17, AAO: 24 },
    { age: 10, ht: 139, wt: 34, SOV: 23, STJ: 18, AAO: 28 }
  ];
  const STORE = "aortic-visits-v1";
  let visits = null;
  try { const s = JSON.parse(localStorage.getItem(STORE)); if (Array.isArray(s) && s.length) visits = s; } catch (e) { /* storage unavailable */ }
  let isExample = !visits;
  if (!visits) visits = EXAMPLE.map((v) => Object.assign({}, v));
  const save = () => { try { localStorage.setItem(STORE, JSON.stringify(visits)); } catch (e) { /* ignore */ } };

  const FIELDS = [["age", "Age (y)"], ["ht", "Height (cm)"], ["wt", "Weight (kg)"], ["SOV", "SOV (mm)"], ["STJ", "STJ (mm)"], ["AAO", "AAo (mm)"]];
  const num = (v) => { const x = parseFloat(String(v == null ? "" : v).replace(",", ".")); return isFinite(x) ? x : NaN; };

  function renderVisitInputs() {
    const tb = $("visit-in"); tb.textContent = "";
    visits.forEach((v, i) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td class="num sub">${i + 1}</td>` + FIELDS.map(([k, lab]) =>
        `<td><input class="cell" id="v${i}-${k}" data-i="${i}" data-k="${k}" inputmode="decimal" autocomplete="off" aria-label="Visit ${i + 1} ${lab}" value="${isFinite(num(v[k])) ? v[k] : ""}"></td>`).join("") +
        `<td><button type="button" class="btn" data-del="${i}" aria-label="Remove visit ${i + 1}">Remove</button></td>`;
      tb.appendChild(tr);
    });
    $("example-note").hidden = !isExample;
  }
  $("visit-in").addEventListener("input", (e) => {
    const t = e.target; if (!t.dataset.k) return;
    visits[+t.dataset.i][t.dataset.k] = t.value; isExample = false; $("example-note").hidden = true; save(); refresh();
  });
  $("visit-in").addEventListener("click", (e) => {
    const i = e.target.dataset.del; if (i == null) return;
    visits.splice(+i, 1); isExample = false; save(); renderVisitInputs(); refresh();
  });
  $("add-visit").addEventListener("click", () => { visits.push({}); isExample = false; save(); renderVisitInputs(); refresh(); const f = $(`v${visits.length - 1}-age`); if (f) f.focus(); });
  $("clear-visits").addEventListener("click", () => { visits = [{}]; isExample = false; save(); renderVisitInputs(); refresh(); });
  $("load-example").addEventListener("click", () => { visits = EXAMPLE.map((v) => Object.assign({}, v)); isExample = true; try { localStorage.removeItem(STORE); } catch (e) { /* ignore */ } renderVisitInputs(); refresh(); });

  function visitResults() {
    return visits.map((v, i) => {
      const p = { ht: num(v.ht), wt: num(v.wt), age: num(v.age), male: state.male };
      const ok = p.ht > 0 && p.age >= 0;
      const hasWt = p.wt > 0;
      const out = { i, p, ok, hasWt, levels: {} };
      if (ok) KEYS.forEach((k) => {
        const d = num(v[k]); if (!(d > 0)) return;
        const r = ref(k, p, hasWt);
        out.levels[k] = { d, z: zOf(d, r), pct: 100 * d / medMm(r) };
      });
      return out;
    });
  }

  function renderVisitResults(res) {
    const tb = $("visit-out"); tb.textContent = "";
    res.forEach((v) => {
      const tr = document.createElement("tr");
      let html = `<td class="num sub">${v.i + 1}</td><td class="num">${v.ok ? v.p.age.toFixed(1) : "—"}</td>`;
      KEYS.forEach((k) => {
        const L = v.levels[k];
        if (!L) { html += `<td class="num sub">—</td><td class="num sub">—</td>`; return; }
        const [cls] = category(L.z);
        html += `<td class="num"><span class="chip ${cls}" style="font-size:12px;padding:1px 7px">${fmtZ(L.z)}</span></td><td class="num">${Math.round(L.pct)}%</td>`;
      });
      html += `<td class="sub">${v.ok ? (v.hasWt ? "Full" : "Height + sex") : "Needs age and height"}</td>`;
      tr.innerHTML = html; tb.appendChild(tr);
    });
  }

  // ---------- chart ----------
  let W = 760, HH = 440;
  const ml = 52, mr = 74, mt = 14, mb = 46;
  const tip = $("tip");
  function curveRef(h) { return ref(state.key, { ht: h, male: state.male }, false); }

  function draw(res) {
    const svg = $("chart"); while (svg.firstChild) svg.removeChild(svg.firstChild);
    const narrow = svg.parentElement.clientWidth < 560;
    W = narrow ? 440 : 760; HH = narrow ? 400 : 440;
    svg.setAttribute("viewBox", `0 0 ${W} ${HH}`);
    const hs = []; for (let h = H0; h <= H1; h += 1) hs.push(h);
    const refs = hs.map(curveRef);
    const lo = Math.min(...refs.map((r) => dAt(r, -3))), hi = Math.max(...refs.map((r) => dAt(r, 6)));
    const ptsD = res.map((v) => v.levels[state.key]).filter(Boolean).map((L) => L.d);
    const y0 = Math.max(0, Math.floor(Math.min(lo, ...ptsD) / 5) * 5), y1 = Math.ceil(Math.max(hi, ...ptsD) / 5) * 5;
    const X = (h) => ml + (h - H0) / (H1 - H0) * (W - ml - mr);
    const Y = (d) => mt + (y1 - Math.min(Math.max(d, y0), y1)) / (y1 - y0) * (HH - mt - mb);

    // bands
    FILLS.forEach(([a, b, c]) => {
      const top = hs.map((h, i) => `${X(h).toFixed(1)},${Y(dAt(refs[i], b)).toFixed(1)}`);
      const bot = hs.map((h, i) => `${X(h).toFixed(1)},${Y(dAt(refs[i], a)).toFixed(1)}`).reverse();
      svg.appendChild(el("polygon", { points: top.concat(bot).join(" "), style: `fill:var(${c});opacity:.75` }));
    });
    // grid + axes
    for (let d = y0; d <= y1; d += 5) {
      svg.appendChild(el("line", { x1: ml, x2: W - mr, y1: Y(d), y2: Y(d), style: "stroke:var(--line);stroke-width:1;opacity:.6" }));
      svg.appendChild(el("text", { x: ml - 8, y: Y(d) + 4, "text-anchor": "end", style: "fill:var(--muted);font:10.5px var(--font-num)" }, d));
    }
    for (let h = 50; h <= 200; h += (W < 600 ? 50 : 25)) {
      if (h < H0 || h > H1) continue;
      svg.appendChild(el("line", { x1: X(h), x2: X(h), y1: HH - mb, y2: HH - mb + 5, style: "stroke:var(--muted)" }));
      svg.appendChild(el("text", { x: X(h), y: HH - mb + 18, "text-anchor": "middle", style: "fill:var(--muted);font:10.5px var(--font-num)" }, h));
    }
    svg.appendChild(el("line", { x1: ml, x2: W - mr, y1: HH - mb, y2: HH - mb, style: "stroke:var(--muted)" }));
    svg.appendChild(el("text", { x: (ml + W - mr) / 2, y: HH - 6, "text-anchor": "middle", style: "fill:var(--muted);font:11.5px var(--font-ui)" }, "Height (cm)"));
    svg.appendChild(el("text", { x: 14, y: (mt + HH - mb) / 2, transform: `rotate(-90 14 ${(mt + HH - mb) / 2})`, "text-anchor": "middle", style: "fill:var(--muted);font:11.5px var(--font-ui)" }, INFO[state.key].name + " (mm)"));
    // curves + direct labels
    ZS.forEach((z) => {
      const d = hs.map((h, i) => `${i ? "L" : "M"}${X(h).toFixed(1)},${Y(dAt(refs[i], z)).toFixed(1)}`).join("");
      svg.appendChild(el("path", { d, style: `fill:none;stroke:var(--ink);stroke-width:${z === 0 ? 2 : 1.2};opacity:${z === 0 ? .9 : .55}` }));
      const yEnd = Y(dAt(refs[refs.length - 1], z));
      svg.appendChild(el("text", { x: W - mr + 6, y: yEnd + 4, style: "fill:var(--ink);font:500 10.5px var(--font-num)" }, z === 0 ? "median" : "z " + (z > 0 ? "+" : "−") + Math.abs(z)));
    });
    // crosshair layer
    const cross = el("line", { y1: mt, y2: HH - mb, style: "stroke:var(--ink);stroke-width:1;opacity:0", "pointer-events": "none" });
    svg.appendChild(cross);
    const hit = el("rect", { x: ml, y: mt, width: W - ml - mr, height: HH - mt - mb, style: "fill:transparent" });
    svg.appendChild(hit);
    hit.addEventListener("mousemove", (e) => {
      const bb = svg.getBoundingClientRect(), sx = (e.clientX - bb.left) * W / bb.width;
      const h = Math.round(H0 + (sx - ml) / (W - ml - mr) * (H1 - H0));
      if (h < H0 || h > H1) return;
      const r = curveRef(h);
      cross.setAttribute("x1", X(h)); cross.setAttribute("x2", X(h)); cross.style.opacity = .35;
      tip.innerHTML = `Height <b>${h} cm</b><br>Median <b>${medMm(r).toFixed(1)} mm</b><br>Normal <b>${dAt(r, -2).toFixed(1)}–${dAt(r, 2).toFixed(1)} mm</b><br>z +3.5 <b>${dAt(r, 3.5).toFixed(1)}</b> · z +5 <b>${dAt(r, 5).toFixed(1)} mm</b>`;
      showTip(e);
    });
    hit.addEventListener("mouseleave", () => { cross.style.opacity = 0; tip.hidden = true; });

    // visits
    const pts = res.filter((v) => v.levels[state.key]).map((v) => ({ v, L: v.levels[state.key], x: X(v.p.ht), y: Y(v.levels[state.key].d) }));
    if (pts.length > 1) svg.appendChild(el("path", { d: pts.map((q, i) => `${i ? "L" : "M"}${q.x},${q.y}`).join(""), style: `fill:none;stroke:var(${INFO[state.key].color});stroke-width:2` }));
    pts.forEach((q) => {
      svg.appendChild(el("circle", { cx: q.x, cy: q.y, r: 5.5, style: `fill:var(${INFO[state.key].color});stroke:var(--surface);stroke-width:2` }));
      svg.appendChild(el("text", { x: q.x - 8, y: q.y - 9, "text-anchor": "end", style: "fill:var(--ink);font:600 10.5px var(--font-num)" }, isFinite(q.v.p.age) ? q.v.p.age.toFixed(q.v.p.age % 1 ? 1 : 0) + " y" : ""));
      const h = el("circle", { cx: q.x, cy: q.y, r: 12, style: "fill:transparent;cursor:default" });
      h.addEventListener("mousemove", (e) => {
        tip.innerHTML = `Visit ${q.v.i + 1} · age <b>${q.v.p.age}</b> y<br>${INFO[state.key].short} <b>${q.L.d} mm</b> at <b>${q.v.p.ht} cm</b><br>z <b>${fmtZ(q.L.z)}</b> · <b>${Math.round(q.L.pct)}%</b> of expected`;
        showTip(e);
      });
      h.addEventListener("mouseleave", () => { tip.hidden = true; });
      svg.appendChild(h);
    });
    $("chart-title").textContent = `${INFO[state.key].name}, ${state.male ? "boys" : "girls"}`;
  }
  function showTip(e) {
    tip.hidden = false;
    const x = Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8);
    tip.style.left = x + "px"; tip.style.top = (e.clientY + 14) + "px";
  }

  // ---------- table view of the curves ----------
  function renderTable() {
    const tb = $("curve-tab"); tb.textContent = "";
    for (let h = 50; h <= 200; h += 10) {
      if (h < H0 || h > H1) continue;
      const r = curveRef(h);
      const tr = document.createElement("tr");
      tr.innerHTML = `<td class="num">${h}</td>` + ZS.map((z) => `<td class="num">${dAt(r, z).toFixed(1)}</td>`).join("");
      tb.appendChild(tr);
    }
    $("curve-cap").textContent = `${INFO[state.key].name}, ${state.male ? "boys" : "girls"}: diameter (mm) at each z, height-and-sex model.`;
  }

  // ---------- controls ----------
  function setPressed(groupId, attr, val) {
    document.querySelectorAll(`#${groupId} button`).forEach((b) => b.setAttribute("aria-pressed", b.dataset[attr] === String(val) ? "true" : "false"));
  }
  document.querySelectorAll("#level-tabs button").forEach((b) => b.addEventListener("click", () => { state.key = b.dataset.k; setPressed("level-tabs", "k", state.key); refresh(); }));
  document.querySelectorAll("#sex-tabs button").forEach((b) => b.addEventListener("click", () => { state.male = b.dataset.m === "1"; setPressed("sex-tabs", "m", state.male ? 1 : 0); refresh(); }));

  function refresh() {
    const res = visitResults();
    draw(res); renderVisitResults(res); renderTable();
  }
  renderVisitInputs();
  refresh();
  let rt; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => draw(visitResults()), 150); });
})();
