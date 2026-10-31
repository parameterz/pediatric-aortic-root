/* Diagnostics page: draws the aggregated summaries in models/diagnostics.js (window.AORTIC_DIAG),
 * written by R/export_diagnostics.R. No individual records are used. */
(function () {
  "use strict";
  const { INFO, el } = AorticModel;
  const $ = (id) => document.getElementById(id);
  const D = window.AORTIC_DIAG;
  if (!D) { $("missing").hidden = false; $("content").hidden = true; return; }
  $("diag-date").textContent = D.meta.created;
  const KEYS = ["SOV", "STJ", "AAO"].filter((k) => D.levels[k]);
  let key = "AAO";

  const tip = $("tip");
  function showTip(e, html) {
    tip.innerHTML = html; tip.hidden = false;
    tip.style.left = Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8) + "px";
    tip.style.top = (e.clientY + 14) + "px";
  }
  const hideTip = () => { tip.hidden = true; };
  const f2 = (v) => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(2);
  const niceTicks = (a, b, n) => {
    const span = b - a, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= n) || 10 * mag;
    const t = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) t.push(+v.toFixed(6));
    return t;
  };

  // Generic panel with linear axes
  function panel(svg, o) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = o.W || 330, H = o.H || 230, m = Object.assign({ l: 38, r: 10, t: 12, b: 38 }, o.m || {});
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
  const line = (pts, style) => el("path", { d: pts.map(([x, y], i) => `${i ? "L" : "M"}${x.toFixed(1)},${y.toFixed(1)}`).join(""), style });
  const area = (top, bot, style) => el("path", { d: top.concat(bot.slice().reverse()).map(([x, y], i) => `${i ? "L" : "M"}${x.toFixed(1)},${y.toFixed(1)}`).join("") + "Z", style });
  const color = () => `var(${INFO[key].color})`;

  // ---------- summary ----------
  function renderSummary(L) {
    $("s-n").textContent = L.n.toLocaleString("en-US");
    $("s-lam").textContent = L.lambda.toFixed(2) + " ± " + L.lambda_se.toFixed(2);
    $("s-skew").textContent = f2(L.moments.skew);
    $("s-kurt").textContent = f2(L.moments.ex_kurt);
    $("s-tails").textContent = L.tails.above[0].toFixed(1) + "% / " + L.tails.below[0].toFixed(1) + "%";
    $("s-sw").textContent = "W " + L.moments.shapiro_W.toFixed(3);
  }

  // ---------- residuals vs covariates ----------
  const COVS = [["height", "Height (cm)"], ["age", "Age (years)"], ["weight", "Weight (kg)"], ["bmi", "BMI (kg/m²)"], ["wfh", "Weight for height (% vs typical)"]];
  function renderResiduals(L) {
    const host = $("resid"); host.textContent = "";
    COVS.forEach(([c, lab]) => {
      const b = L.bins[c], t = L.trends[c];
      const fig = document.createElement("figure"); fig.className = "mini";
      const svg = el("svg", { role: "img", "aria-label": `z-score against ${lab}` });
      fig.appendChild(svg);
      const cap = document.createElement("figcaption"); cap.textContent = lab; fig.appendChild(cap);
      host.appendChild(fig);
      const x0 = Math.min(...b.xlo), x1 = Math.max(...b.xhi);
      const P = panel(svg, { x0, x1, y0: -3, y1: 3, yt: [-3, -2, -1, 0, 1, 2, 3], ylab: "z", yfmt: (v) => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) });
      const xs = b.x.map(P.X);
      svg.appendChild(area(xs.map((x, i) => [x, P.Y(b.q90[i])]), xs.map((x, i) => [x, P.Y(b.q10[i])]), `fill:${color()};opacity:.14`));
      svg.appendChild(area(xs.map((x, i) => [x, P.Y(b.q75[i])]), xs.map((x, i) => [x, P.Y(b.q25[i])]), `fill:${color()};opacity:.22`));
      [-2, 2].forEach((z) => svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(z), y2: P.Y(z), style: "stroke:var(--muted);stroke-dasharray:3 3" })));
      svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(0), y2: P.Y(0), style: "stroke:var(--muted)" }));
      const tx = t.x.map(P.X);
      svg.appendChild(area(tx.map((x, i) => [x, P.Y(t.fit[i] + 2 * t.se[i])]), tx.map((x, i) => [x, P.Y(t.fit[i] - 2 * t.se[i])]), "fill:var(--ink);opacity:.12"));
      svg.appendChild(line(tx.map((x, i) => [x, P.Y(t.fit[i])]), "fill:none;stroke:var(--ink);stroke-width:2"));
      b.x.forEach((xv, i) => {
        svg.appendChild(el("circle", { cx: P.X(xv), cy: P.Y(b.mean[i]), r: 2.8, style: "fill:var(--ink)" }));
        const hit = el("rect", { x: P.X(b.xlo[i]), y: P.m.t, width: Math.max(4, P.X(b.xhi[i]) - P.X(b.xlo[i])), height: P.H - P.m.t - P.m.b, style: "fill:transparent" });
        hit.addEventListener("mousemove", (e) => showTip(e,
          `${lab}: <b>${b.xlo[i].toFixed(1)}–${b.xhi[i].toFixed(1)}</b><br>n <b>${b.n[i]}</b> · mean z <b>${f2(b.mean[i])}</b> · SD <b>${b.sd[i].toFixed(2)}</b><br>above +2 <b>${b.gt2[i].toFixed(1)}%</b> · below −2 <b>${b.lt2[i].toFixed(1)}%</b>`));
        hit.addEventListener("mouseleave", hideTip);
        svg.appendChild(hit);
      });
    });
    // sex
    const s = L.sex, tb = $("sex-tab"); tb.textContent = "";
    s.sex.forEach((sx, i) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td>${sx === "M" ? "Male" : "Female"}</td><td class="num">${s.n[i]}</td><td class="num">${f2(s.mean[i])}</td><td class="num">${(2 * s.se[i]).toFixed(2)}</td><td class="num">${s.sd[i].toFixed(2)}</td><td class="num">${s.gt2[i].toFixed(1)}%</td><td class="num">${s.lt2[i].toFixed(1)}%</td>`;
      tb.appendChild(tr);
    });
  }

  // ---------- Q-tests ----------
  function renderQ(L) {
    const q = L.qstats, tb = $("q-tab"); tb.textContent = "";
    const cols = ["Z1", "Z2", "Z3", "Z4"];
    q.group.forEach((g, i) => {
      const isTot = /TOTAL|df for|p-val/.test(g);
      const tr = document.createElement("tr");
      if (/p-val/.test(g)) tr.className = "primary";
      let html = `<td>${g.replace("TOTAL Q stats", "Total Q").replace("df for Q stats", "Degrees of freedom").replace("p-val for Q stats", "p-value")}</td>`;
      cols.forEach((c) => {
        const v = q[c][i];
        let cls = "";
        if (!isTot && Math.abs(v) > 2) cls = "chip mild";
        if (/p-val/.test(g) && v < 0.05) cls = "chip mild";
        const txt = /p-val/.test(g) ? v.toFixed(3) : /df for/.test(g) ? v.toFixed(1) : isTot ? v.toFixed(1) : f2(v);
        html += `<td class="num">${cls ? `<span class="${cls}" style="font-size:12px;padding:1px 7px">${txt}</span>` : txt}</td>`;
      });
      html += `<td class="num sub">${isTot ? "" : q.N[i]}</td>`;
      tr.innerHTML = html; tb.appendChild(tr);
    });
  }

  // ---------- distribution ----------
  function renderDist(L) {
    // Q-Q
    const q = L.qq;
    const lim = Math.ceil(Math.max(3, ...q.samp.map(Math.abs), ...q.theo.map(Math.abs)) * 2) / 2;
    let P = panel($("qq"), { W: 330, H: 300, x0: -lim, x1: lim, y0: -lim, y1: lim, xlab: "Expected z (normal)", ylab: "Observed z" });
    $("qq").appendChild(line([[P.X(-lim), P.Y(-lim)], [P.X(lim), P.Y(lim)]], "stroke:var(--muted);stroke-dasharray:4 3"));
    q.p.forEach((p, i) => {
      $("qq").appendChild(el("circle", { cx: P.X(q.theo[i]), cy: P.Y(q.samp[i]), r: 3.5, style: `fill:${color()};stroke:var(--surface);stroke-width:1.5` }));
      const hit = el("circle", { cx: P.X(q.theo[i]), cy: P.Y(q.samp[i]), r: 8, style: "fill:transparent" });
      hit.addEventListener("mousemove", (e) => showTip(e, `Percentile <b>${(100 * p).toFixed(1)}</b><br>expected <b>${f2(q.theo[i])}</b> · observed <b>${f2(q.samp[i])}</b>`));
      hit.addEventListener("mouseleave", hideTip);
      $("qq").appendChild(hit);
    });
    // histogram + normal curve
    const h = L.hist, n = h.counts.reduce((a, b) => a + b, 0), bw = h.breaks[1] - h.breaks[0];
    const x0 = Math.min(-4, h.breaks[0]), x1 = Math.max(4, h.breaks[h.breaks.length - 1]);
    const ymax = Math.max(...h.counts, n * bw / Math.sqrt(2 * Math.PI)) * 1.08;
    P = panel($("hist"), { W: 330, H: 300, x0, x1, y0: 0, y1: ymax, xlab: "z", ylab: "Children" });
    h.counts.forEach((c, i) => {
      const xa = P.X(h.breaks[i]) + 1, xb = P.X(h.breaks[i + 1]) - 1;
      const r = el("rect", { x: xa, y: P.Y(c), width: Math.max(1, xb - xa), height: P.Y(0) - P.Y(c), rx: 1.5, style: `fill:${color()};opacity:.75` });
      r.addEventListener("mousemove", (e) => showTip(e, `z <b>${f2(h.breaks[i])}</b> to <b>${f2(h.breaks[i + 1])}</b><br><b>${c}</b> children (normal: ${(n * (Phi(h.breaks[i + 1]) - Phi(h.breaks[i]))).toFixed(1)})`));
      r.addEventListener("mouseleave", hideTip);
      $("hist").appendChild(r);
    });
    const xs = []; for (let v = x0; v <= x1 + 1e-9; v += 0.05) xs.push(v);
    $("hist").appendChild(line(xs.map((v) => [P.X(v), P.Y(n * bw * Math.exp(-v * v / 2) / Math.sqrt(2 * Math.PI))]), "fill:none;stroke:var(--ink);stroke-width:2"));
    // tails
    const t = L.tails, tb = $("tail-tab"); tb.textContent = "";
    t.cut.forEach((c, i) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td class="num">±${c}</td><td class="num">${t.expected[i].toFixed(2)}%</td><td class="num">${t.above[i].toFixed(2)}%</td><td class="num">${t.below[i].toFixed(2)}%</td><td class="num">${(t.expected[i] / 100 * L.n).toFixed(0)} / ${(t.above[i] / 100 * L.n).toFixed(0)} / ${(t.below[i] / 100 * L.n).toFixed(0)}</td>`;
      tb.appendChild(tr);
    });
  }
  // normal CDF (Abramowitz-Stegun 7.1.26)
  function Phi(x) {
    const t = 1 / (1 + 0.3275911 * Math.abs(x) / Math.SQRT2);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x / 2);
    return x >= 0 ? (1 + y) / 2 : (1 - y) / 2;
  }

  // ---------- worm plots ----------
  function renderWorm(L) {
    const host = $("worm"); host.textContent = "";
    L.worm.forEach((g) => {
      const fig = document.createElement("figure"); fig.className = "mini";
      const svg = el("svg", { role: "img", "aria-label": `Worm plot, ${g.label}` });
      fig.appendChild(svg);
      const cap = document.createElement("figcaption"); cap.textContent = `${g.label} (n ${g.n})`; fig.appendChild(cap);
      host.appendChild(fig);
      const P = panel(svg, { W: 330, H: 200, x0: -2.5, x1: 2.5, y0: -0.6, y1: 0.6, yt: [-0.6, -0.3, 0, 0.3, 0.6], xlab: "Expected z", ylab: "Deviation" });
      const band = g.theo.map((q) => { const p = Phi(q), dens = Math.exp(-q * q / 2) / Math.sqrt(2 * Math.PI); return 1.96 * Math.sqrt(p * (1 - p) / g.n) / dens; });
      const xs = g.theo.map(P.X);
      svg.appendChild(area(xs.map((x, i) => [x, P.Y(band[i])]), xs.map((x, i) => [x, P.Y(-band[i])]), "fill:var(--ink);opacity:.08"));
      svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(0), y2: P.Y(0), style: "stroke:var(--muted)" }));
      svg.appendChild(line(xs.map((x, i) => [x, P.Y(g.dev[i])]), `fill:none;stroke:${color()};stroke-width:2`));
      g.dev.forEach((dv, i) => svg.appendChild(el("circle", { cx: xs[i], cy: P.Y(dv), r: 2.2, style: `fill:${Math.abs(dv) > band[i] ? "var(--chip-severe-ink)" : color()}` })));
    });
  }

  // ---------- shape comparison (all levels) ----------
  function renderShape() {
    const tb = $("aic-tab"); tb.textContent = "";
    KEYS.forEach((k) => {
      const L = D.levels[k], a = L.aic, best = a.free;
      const cell = (v) => v == null ? "—" : (v - best >= 0 ? "+" : "−") + Math.abs(v - best).toFixed(1);
      const tr = document.createElement("tr");
      tr.innerHTML = `<td style="white-space:nowrap"><span class="dot" style="background:var(${INFO[k].color})"></span>${INFO[k].name}</td><td class="num">${L.lambda.toFixed(2)} ± ${L.lambda_se.toFixed(2)}</td>` +
        `<td class="num">0</td><td class="num">${cell(a.log)}</td><td class="num">${cell(a.linear)}</td><td class="num">${cell(a.squared)}</td><td class="num">${cell(a.bct)}</td>`;
      tb.appendChild(tr);
    });
  }

  // ---------- local exponent (all levels) ----------
  function renderExponent() {
    const svg = $("expo");
    const all = KEYS.map((k) => D.levels[k].exponent);
    const x0 = Math.min(...all.map((e) => e.ht[0])), x1 = Math.max(...all.map((e) => e.ht[e.ht.length - 1]));
    const P = panel(svg, { W: 700, H: 330, m: { l: 44, r: 130, t: 14, b: 42 }, x0, x1, y0: 0.4, y1: 1.3, yt: [0.5, 0.75, 1, 1.25], xlab: "Height (cm)", ylab: "Local exponent" });
    [[1, "geometric (1)"], [0.5, "adult-like (0.5)"]].forEach(([v, lab]) => {
      svg.appendChild(el("line", { x1: P.m.l, x2: P.W - P.m.r, y1: P.Y(v), y2: P.Y(v), style: "stroke:var(--muted);stroke-dasharray:4 3" }));
      svg.appendChild(el("text", { x: P.m.l + 4, y: P.Y(v) - 4, style: "fill:var(--muted);font:10px var(--font-num)" }, lab));
    });
    const ends = [];
    KEYS.forEach((k) => {
      const e = D.levels[k].exponent, xs = e.ht.map(P.X);
      svg.appendChild(area(xs.map((x, i) => [x, P.Y(e.b[i] + 2 * e.se[i])]), xs.map((x, i) => [x, P.Y(e.b[i] - 2 * e.se[i])]), `fill:var(${INFO[k].color});opacity:.12`));
      svg.appendChild(line(xs.map((x, i) => [x, P.Y(e.b[i])]), `fill:none;stroke:var(${INFO[k].color});stroke-width:2`));
      ends.push({ y: P.Y(e.b[e.b.length - 1]), k });
    });
    ends.sort((a, b) => a.y - b.y);
    for (let i = 1; i < ends.length; i++) if (ends[i].y - ends[i - 1].y < 13) ends[i].y = ends[i - 1].y + 13;
    ends.forEach((e) => svg.appendChild(el("text", { x: P.W - P.m.r + 8, y: e.y + 4, style: "fill:var(--ink);font:11px var(--font-ui)" }, INFO[e.k].name)));
    const hit = el("rect", { x: P.m.l, y: P.m.t, width: P.W - P.m.l - P.m.r, height: P.H - P.m.t - P.m.b, style: "fill:transparent" });
    hit.addEventListener("mousemove", (ev) => {
      const bb = svg.getBoundingClientRect(), sx = (ev.clientX - bb.left) * P.W / bb.width;
      const h = x0 + (sx - P.m.l) / (P.W - P.m.l - P.m.r) * (x1 - x0);
      showTip(ev, `Height <b>${h.toFixed(0)} cm</b><br>` + KEYS.map((k) => {
        const e = D.levels[k].exponent; let i = e.ht.findIndex((v) => v >= h); if (i < 0) i = e.ht.length - 1;
        return `${INFO[k].short} <b>${e.b[i].toFixed(2)}</b> ± ${(2 * e.se[i]).toFixed(2)}`;
      }).join("<br>"));
    });
    hit.addEventListener("mouseleave", hideTip);
    svg.appendChild(hit);
  }

  function renderLevel() {
    const L = D.levels[key];
    document.querySelectorAll(".lvl-name").forEach((s) => (s.textContent = INFO[key].name));
    renderSummary(L); renderResiduals(L); renderQ(L); renderDist(L); renderWorm(L);
  }
  document.querySelectorAll("#level-tabs button").forEach((b) => b.addEventListener("click", () => {
    key = b.dataset.k;
    document.querySelectorAll("#level-tabs button").forEach((x) => x.setAttribute("aria-pressed", x === b ? "true" : "false"));
    renderLevel();
  }));
  renderLevel(); renderShape(); renderExponent();
})();
