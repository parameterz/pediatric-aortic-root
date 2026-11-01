/* Method comparison page: draws the aggregated summaries in models/methods.js (window.AORTIC_METHODS),
 * written by R/export_methods_comparison.R. No individual records are used. */
(function () {
  "use strict";
  const { INFO, el } = AorticModel;
  const $ = (id) => document.getElementById(id);
  const D = window.AORTIC_METHODS;
  if (!D) { $("missing").hidden = false; $("content").hidden = true; return; }
  $("m-date").textContent = D.meta.created;

  // Fixed color per method (categorical slots in order; never re-assigned when filtering)
  const COLOR = { ours_cv: "--m1", phn: "--m2", ls: "--m3", csa_tab: "--m4", mah_ht: "--m5", hs: "--m6", mah_bsa: "--m7", csa_sheet: "--m8", csa_refit: "--m8" };
  const DASHED = new Set(["csa_refit"]);   // same family as csa_sheet: same hue, dashed (no ninth color)
  const NAME = {
    ours_cv: "This model, cross-validated", ours: "This model, in sample", hs: "Height and sex (this model, no weight)",
    phn: "PHN z-score (BSA⁰·⁵)", mah_bsa: "Mahgerefteh BSA⁰·⁵ (≥ 2 y)", mah_ht: "Mahgerefteh height (≥ 2 y)",
    ls: "Lauzon-Schnittka", csa_tab: "CSA/ht, age-group table", csa_sheet: "CSA/ht, continuous (as distributed)",
    csa_refit: "CSA/ht, continuous (refit)"
  };
  const SHORT = { ours_cv: "This model (CV)", hs: "Height and sex", phn: "PHN BSA", mah_bsa: "Mahgerefteh BSA", mah_ht: "Mahgerefteh height", ls: "Lauzon-Schnittka", csa_tab: "CSA/ht table", csa_sheet: "CSA/ht as distributed", csa_refit: "CSA/ht refit" };
  const col = (m) => COLOR[m] ? `var(${COLOR[m]})` : "var(--muted)";
  const swBg = (m) => DASHED.has(m) ? `repeating-linear-gradient(90deg,${col(m)} 0 4px,transparent 4px 7px)` : col(m);
  // Display order, in groups: this model, its no-weight version, BSA, Mahgerefteh, CSA/ht, Lauzon-Schnittka
  const GROUPS = [["ours_cv", "ours"], ["hs"], ["phn"], ["mah_ht", "mah_bsa"], ["csa_tab", "csa_sheet", "csa_refit"], ["ls"]];
  const ORDER = GROUPS.flat();
  const GROUP_OF = {}; GROUPS.forEach((g, i) => g.forEach((m) => { GROUP_OF[m] = i; }));
  const byOrder = (a, b) => ORDER.indexOf(a) - ORDER.indexOf(b);
  const CHARTABLE = ORDER.filter((m) => COLOR[m]);
  const selected = new Set(["ours_cv", "phn", "ls", "csa_tab"]);
  let key = "SOV", stat = "mean";

  const tip = $("tip");
  function showTip(e, html) {
    tip.innerHTML = html; tip.hidden = false;
    tip.style.left = Math.max(8, Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8)) + "px";
    tip.style.top = (e.clientY + 14) + "px";
  }
  const hideTip = () => { tip.hidden = true; };
  const sgn = (v, d) => { const r = +v.toFixed(d); return (r > 0 ? "+" : r < 0 ? "−" : "") + Math.abs(r).toFixed(d); };
  const niceTicks = (a, b, n) => {
    const span = b - a, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= n) || 10 * mag;
    const t = []; for (let v = Math.ceil(a / step - 1e-9) * step; v <= b + 1e-9; v += step) t.push(+v.toFixed(6));
    return t;
  };
  function panel(svg, o) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = o.W || 330, H = o.H || 230, m = Object.assign({ l: 40, r: 10, t: 12, b: 38 }, o.m || {});
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const X = (v) => m.l + (v - o.x0) / (o.x1 - o.x0) * (W - m.l - m.r);
    const Y = (v) => m.t + (o.y1 - Math.min(Math.max(v, o.y0), o.y1)) / (o.y1 - o.y0) * (H - m.t - m.b);
    (o.yt || niceTicks(o.y0, o.y1, 5)).forEach((v) => {
      svg.appendChild(el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), style: "stroke:var(--line);opacity:.8" }));
      svg.appendChild(el("text", { x: m.l - 6, y: Y(v) + 3.5, "text-anchor": "end", style: "fill:var(--muted);font:10px var(--font-num)" }, o.yfmt ? o.yfmt(v) : v));
    });
    (o.xt || niceTicks(o.x0, o.x1, 5)).forEach((v) => {
      svg.appendChild(el("line", { x1: X(v), x2: X(v), y1: H - m.b, y2: H - m.b + 4, style: "stroke:var(--muted)" }));
      svg.appendChild(el("text", { x: X(v), y: H - m.b + 15, "text-anchor": "middle", style: "fill:var(--muted);font:10px var(--font-num)" }, o.xfmt ? o.xfmt(v) : v));
    });
    svg.appendChild(el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, style: "stroke:var(--muted)" }));
    if (o.xlab) svg.appendChild(el("text", { x: (m.l + W - m.r) / 2, y: H - 5, "text-anchor": "middle", style: "fill:var(--muted);font:11px var(--font-ui)" }, o.xlab));
    if (o.ylab) svg.appendChild(el("text", { x: 11, y: (m.t + H - m.b) / 2, transform: `rotate(-90 11 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", style: "fill:var(--muted);font:11px var(--font-ui)" }, o.ylab));
    return { X, Y, W, H, m };
  }
  // polyline with gaps at null values
  function path(xs, ys, P, style) {
    let d = "", pen = false;
    xs.forEach((x, i) => {
      const y = ys[i];
      if (y === null || y === undefined || !isFinite(y)) { pen = false; return; }
      d += `${pen ? "L" : "M"}${P.X(x).toFixed(1)},${P.Y(y).toFixed(1)}`; pen = true;
    });
    return el("path", { d, style: "fill:none;stroke-linejoin:round;stroke-linecap:round;" + style });
  }
  // invisible hover layer: nearest x index
  function hover(svg, P, xs, html) {
    const r = el("rect", { x: P.m.l, y: P.m.t, width: P.W - P.m.l - P.m.r, height: P.H - P.m.t - P.m.b, style: "fill:transparent" });
    const cross = el("line", { y1: P.m.t, y2: P.H - P.m.b, style: "stroke:var(--muted);stroke-dasharray:2 3;display:none" });
    svg.appendChild(cross); svg.appendChild(r);
    const move = (e) => {
      const pt = svg.getBoundingClientRect(), vx = (e.clientX - pt.left) / pt.width * P.W;
      let best = 0, bd = Infinity;
      xs.forEach((x, i) => { const dd = Math.abs(P.X(x) - vx); if (dd < bd) { bd = dd; best = i; } });
      cross.setAttribute("x1", P.X(xs[best])); cross.setAttribute("x2", P.X(xs[best])); cross.style.display = "";
      showTip(e, html(best));
    };
    r.addEventListener("pointermove", move);
    r.addEventListener("pointerleave", () => { cross.style.display = "none"; hideTip(); });
  }

  // ---------- scorecard ----------
  const T = {
    mean: (v) => grade(Math.abs(v), 0.1, 0.25),
    sd: (v) => grade(Math.abs(v - 1), 0.05, 0.1),
    skew: (v) => grade(Math.abs(v), 0.1, 0.25),
    tail: (v) => (v < 1.2 || v > 4.0) ? "bad" : (v < 1.6 || v > 3.0) ? "warn" : "",
    slope: (v) => grade(Math.abs(v), 0.05, 0.15)
  };
  function grade(a, w, b) { return a > b ? "bad" : a > w ? "warn" : ""; }
  function renderScore() {
    const tb = $("score-body"); tb.textContent = "";
    D.summary[key].slice().sort((a, b) => byOrder(a.method, b.method)).forEach((r) => {
      const tr = document.createElement("tr");
      if (r.method === "ours_cv" || r.method === "ours") tr.className = "ours";
      const cell = (txt, cls) => `<td class="num${cls ? " " + cls : ""}">${txt}</td>`;
      const sw = COLOR[r.method] ? `<span class="sw" style="background:${swBg(r.method)}"></span>` : `<span class="sw"></span>`;
      tr.innerHTML = `<td>${sw}${NAME[r.method] || r.method}</td>` + cell(r.n.toLocaleString("en-US")) +
        cell(sgn(r.mean, 2), T.mean(r.mean)) + cell(r.sd.toFixed(2), T.sd(r.sd)) + cell(sgn(r.skew, 2), T.skew(r.skew)) +
        cell(r.gt2.toFixed(1), T.tail(r.gt2)) + cell(r.lt2.toFixed(1), T.tail(r.lt2)) +
        ["sl_ht", "sl_age", "sl_wfh", "sl_bmiz", "sex_diff"].map((k) => cell(sgn(r[k], 2), T.slope(r[k]))).join("");
      tb.appendChild(tr);
    });
  }

  // ---------- trend panels ----------
  function renderPills() {
    const host = $("method-pills"); host.textContent = "";
    const have = new Set(D.summary[key].map((r) => r.method));
    let grp = null, prev = null;
    CHARTABLE.filter((m) => have.has(m)).forEach((m) => {
      if (GROUP_OF[m] !== prev) {                       // each group stays together when the row wraps
        grp = document.createElement("span"); grp.className = "pill-group"; host.appendChild(grp); prev = GROUP_OF[m];
      }
      const b = document.createElement("button");
      b.type = "button"; b.setAttribute("aria-pressed", selected.has(m) ? "true" : "false");
      b.innerHTML = `<span class="sw" style="background:${swBg(m)}"></span>${SHORT[m]}`;
      b.addEventListener("click", () => { selected.has(m) ? selected.delete(m) : selected.add(m); renderPills(); renderTrends(); });
      grp.appendChild(b);
    });
  }
  const XF = { ht_cm: (v) => v.toFixed(0), age: (v) => v, wfh: (v) => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v), bmiz: (v) => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) };
  const XT = { ht_cm: (v) => v.toFixed(0) + " cm", age: (v) => v.toFixed(1) + " y", wfh: (v) => sgn(v, 0) + "%", bmiz: (v) => sgn(v, 1) };
  function renderTrends() {
    const host = $("trends"); host.textContent = "";
    const have = new Set(D.summary[key].map((r) => r.method));
    const ms = CHARTABLE.filter((m) => have.has(m) && selected.has(m));
    D.bins[key].forEach((b) => {
      const fig = document.createElement("figure"); fig.className = "mini";
      const svg = el("svg", { role: "img", "aria-label": `${stat === "mean" ? "Mean" : "SD of"} z against ${b.label}` });
      fig.appendChild(svg);
      const cap = document.createElement("figcaption"); cap.textContent = b.label; fig.appendChild(cap);
      host.appendChild(fig);
      const S = b[stat];
      let lo = stat === "mean" ? -0.6 : 0.7, hi = stat === "mean" ? 0.6 : 1.3;
      ms.forEach((m) => S[m].forEach((v) => { if (v !== null && isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }));
      const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
      const x0 = Math.min(...b.xlo), x1 = Math.max(...b.xhi);
      const P = panel(svg, { x0, x1, y0: lo, y1: hi, ylab: stat === "mean" ? "Mean z" : "SD of z", xfmt: XF[b.var],
        yfmt: (v) => stat === "mean" ? (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) : v });
      const target = stat === "mean" ? 0 : 1;
      svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(target), y2: P.Y(target), style: "stroke:var(--muted);stroke-width:1.2" }));
      // draw the cross-validated model last so it sits on top
      ms.slice().sort((a, c) => (a === "ours_cv") - (c === "ours_cv")).forEach((m) => {
        svg.appendChild(path(b.x, S[m], P, `stroke:var(--surface);stroke-width:4.5`));
        svg.appendChild(path(b.x, S[m], P, `stroke:${col(m)};stroke-width:${m === "ours_cv" ? 2.5 : 2}${DASHED.has(m) ? ";stroke-dasharray:5 3" : ""}`));
      });
      hover(svg, P, b.x, (i) => `<div style="margin-bottom:4px">${XT[b.var](b.xlo[i])} to ${XT[b.var](b.xhi[i])} · ${b.n[i]} children</div>` +
        ms.map((m) => `<div><span class="sw" style="display:inline-block;width:10px;height:3px;border-radius:2px;margin-right:5px;vertical-align:middle;background:${swBg(m)}"></span>${SHORT[m]}: <b>${S[m][i] === null ? "n/a" : stat === "mean" ? sgn(S[m][i], 2) : S[m][i].toFixed(2)}</b></div>`).join(""));
    });
  }

  // ---------- agreement ----------
  function renderAgree() {
    const tb = $("agree-body"); tb.textContent = "";
    const A = D.agree[key].filter((r) => r.method !== "ours").sort((a, b) => byOrder(a.method, b.method));
    const maxF = Math.max(...A.map((r) => r.both + r.ours_only + r.other_only));
    A.forEach((r) => {
      const tr = document.createElement("tr");
      const seg = (v, c) => v > 0 ? `<span style="width:${(100 * v / maxF).toFixed(2)}%;background:${c}" title="${v}"></span>` : "";
      tr.innerHTML = `<td>${NAME[r.method] || r.method}</td><td class="num">${r.n.toLocaleString("en-US")}</td><td class="num">${r.both}</td><td class="num">${r.ours_only}</td><td class="num">${r.other_only}</td>` +
        `<td><div class="agbar" role="img" aria-label="${r.both} both, ${r.ours_only} this model only, ${r.other_only} other only">${seg(r.ours_only, "var(--m1)")}${seg(r.both, "var(--muted)")}${seg(r.other_only, "var(--m2)")}</div></td>`;
      tb.appendChild(tr);
    });
    if (!$("ag-key")) {
      const k = document.createElement("div"); k.className = "ag-key"; k.id = "ag-key";
      k.innerHTML = `<span><i style="background:var(--m1)"></i>This model only</span><span><i style="background:var(--muted)"></i>Both</span><span><i style="background:var(--m2)"></i>Other method only</span>`;
      $("agree-body").closest(".table-wrap").after(k);
    }
  }

  // ---------- weight curves ----------
  const CURVE = [["ours", "This model", "--m1"], ["ls", "Lauzon-Schnittka", "--m3"], ["phn", "PHN BSA⁰·⁵", "--m2"]];
  function renderCurves() {
    $("curve-legend").innerHTML = CURVE.map(([, n, c]) => `<span class="lp"><span class="sw" style="background:var(${c})"></span>${n}</span>`).join("") +
      `<span class="lp"><span class="sw" style="background:repeating-linear-gradient(90deg,var(--m1) 0 4px,transparent 4px 7px)"></span>This model, beyond the PHN build range (held flat)</span>`;
    const host = $("curves"); host.textContent = "";
    let lo = 95, hi = 105;
    D.curves[key].forEach((c) => CURVE.forEach(([m]) => c.pct[m].forEach((v) => { lo = Math.min(lo, v); hi = Math.max(hi, v); })));
    lo = Math.floor(lo / 5) * 5; hi = Math.ceil(hi / 5) * 5;
    D.curves[key].forEach((c) => {
      const fig = document.createElement("figure"); fig.className = "mini";
      const svg = el("svg", { role: "img", "aria-label": `Expected diameter against BMI-for-age z: ${c.label}` });
      fig.appendChild(svg);
      const cap = document.createElement("figcaption"); cap.textContent = c.label; fig.appendChild(cap);
      host.appendChild(fig);
      const P = panel(svg, { x0: -2, x1: 4, y0: lo, y1: hi, xt: [-2, -1, 0, 1, 2, 3, 4], xlab: "BMI-for-age z (WHO)", ylab: "% of value at BMIz 0",
        xfmt: (v) => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v), yfmt: (v) => v + "%" });
      svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(100), y2: P.Y(100), style: "stroke:var(--muted);stroke-width:1.2" }));
      CURVE.slice().reverse().forEach(([m, , cv]) => {
        if (m === "ours") {
          const inR = c.pct.ours.map((v, i) => (c.beyond_range[i] ? null : v));
          const outR = c.pct.ours.map((v, i) => (c.beyond_range[i] || (c.beyond_range[i - 1] && !c.beyond_range[i]) || c.beyond_range[i + 1] ? v : null));
          svg.appendChild(path(c.bmiz, outR, P, `stroke:var(${cv});stroke-width:2;stroke-dasharray:4 3`));
          svg.appendChild(path(c.bmiz, inR, P, `stroke:var(--surface);stroke-width:4.5`));
          svg.appendChild(path(c.bmiz, inR, P, `stroke:var(${cv});stroke-width:2.5`));
        } else {
          svg.appendChild(path(c.bmiz, c.pct[m], P, `stroke:var(--surface);stroke-width:4.5`));
          svg.appendChild(path(c.bmiz, c.pct[m], P, `stroke:var(${cv});stroke-width:2`));
        }
      });
      hover(svg, P, c.bmiz, (i) => `<div style="margin-bottom:4px">BMIz ${sgn(c.bmiz[i], 2)} · ${c.wt[i].toFixed(1)} kg (${sgn(c.wfh_pct[i], 0)}% vs typical for height)</div>` +
        CURVE.map(([m, n, cv]) => `<div><span style="display:inline-block;width:10px;height:3px;border-radius:2px;margin-right:5px;vertical-align:middle;background:var(${cv})"></span>${n}: <b>${c.mm[m][i].toFixed(1)} mm</b> (${c.pct[m][i].toFixed(1)}%)${m === "ours" && c.beyond_range[i] ? " held flat" : ""}</div>`).join(""));
    });
  }

  function render() {
    document.querySelectorAll(".lvl-name").forEach((s) => { s.textContent = INFO[key].name; });
    renderScore(); renderPills(); renderTrends(); renderAgree(); renderCurves();
  }
  document.querySelectorAll("#level-tabs button").forEach((b) => b.addEventListener("click", () => {
    key = b.dataset.k;
    document.querySelectorAll("#level-tabs button").forEach((x) => x.setAttribute("aria-pressed", x === b ? "true" : "false"));
    render();
  }));
  document.querySelectorAll("#stat-tabs button").forEach((b) => b.addEventListener("click", () => {
    stat = b.dataset.s;
    document.querySelectorAll("#stat-tabs button").forEach((x) => x.setAttribute("aria-pressed", x === b ? "true" : "false"));
    renderTrends();
  }));
  render();
})();
