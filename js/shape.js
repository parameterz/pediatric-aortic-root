/* Distribution shape (lambda) explainer.
 * Holds the fitted median M and spread sigma for the chosen child and level, and lets the
 * Box-Cox power lambda vary, to show how the tail and the z-score of one diameter respond. */
(function () {
  "use strict";
  const { M, KEYS, INFO, ref, zOf, dAt, medMm, category, fmtZ, el } = AorticModel;
  const $ = (id) => document.getElementById(id);
  const state = { key: "AAO", lam: null };
  const num = (id) => { const v = parseFloat(String($(id).value).replace(",", ".")); return isFinite(v) ? v : NaN; };

  function child() {
    const p = { ht: num("ht"), wt: num("wt"), age: num("age"), male: $("sex-m").checked };
    const ok = p.ht > 0 && p.age >= 0;
    return { p, ok, hasWt: p.wt > 0 };
  }
  const fitted = () => M.structures[state.key].full.nu;
  // BCCG density on the mm scale (no truncation correction; sigma is small)
  function density(dMm, r, nu) {
    const z = zOf(dMm, r, nu), y = dMm / 10, mu = Math.exp(r.lm);
    const phi = Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
    return Math.pow(y, nu - 1) / (Math.pow(mu, nu) * r.sigma) * phi / 10;   // per mm
  }

  const tip = $("tip");
  function showTip(e, html) {
    tip.innerHTML = html; tip.hidden = false;
    tip.style.left = Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8) + "px";
    tip.style.top = (e.clientY + 14) + "px";
  }
  function axes(svg, W, H, m, X, Y, xt, yt, xlab, ylab) {
    yt.forEach((v) => {
      svg.appendChild(el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), style: "stroke:var(--line);opacity:.7" }));
      svg.appendChild(el("text", { x: m.l - 8, y: Y(v) + 4, "text-anchor": "end", style: "fill:var(--muted);font:10.5px var(--font-num)" }, v));
    });
    xt.forEach((v) => {
      svg.appendChild(el("line", { x1: X(v), x2: X(v), y1: H - m.b, y2: H - m.b + 5, style: "stroke:var(--muted)" }));
      svg.appendChild(el("text", { x: X(v), y: H - m.b + 18, "text-anchor": "middle", style: "fill:var(--muted);font:10.5px var(--font-num)" }, v));
    });
    svg.appendChild(el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, style: "stroke:var(--muted)" }));
    svg.appendChild(el("text", { x: (m.l + W - m.r) / 2, y: H - 6, "text-anchor": "middle", style: "fill:var(--muted);font:11.5px var(--font-ui)" }, xlab));
    if (ylab) svg.appendChild(el("text", { x: 13, y: (m.t + H - m.b) / 2, transform: `rotate(-90 13 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", style: "fill:var(--muted);font:11.5px var(--font-ui)" }, ylab));
  }
  const ticks = (a, b, step) => { const t = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) t.push(+v.toFixed(2)); return t; };
  const LAMS = () => [[0, "λ 0 (log)"], [fitted(), `λ ${fitted().toFixed(2)} (fitted)`], [1, "λ 1 (linear)"], [2, "λ 2 (squared)"]];

  function drawDensity(r, D) {
    const svg = $("dens"); while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = 640, H = 280, m = { l: 20, r: 20, t: 22, b: 44 };
    const x0 = Math.floor(dAt(r, -4, 0)), x1 = Math.ceil(dAt(r, 6, 0));
    const xs = []; for (let i = 0; i <= 300; i++) xs.push(x0 + (x1 - x0) * i / 300);
    const dSel = xs.map((x) => density(x, r, state.lam)), dFit = xs.map((x) => density(x, r, fitted()));
    const ymax = Math.max(...dSel, ...dFit) * 1.08;
    const X = (v) => m.l + (Math.min(Math.max(v, x0), x1) - x0) / (x1 - x0) * (W - m.l - m.r), Y = (v) => m.t + (1 - v / ymax) * (H - m.t - m.b);
    axes(svg, W, H, m, X, Y, ticks(x0, x1, x1 - x0 > 30 ? 10 : 5), [], "Diameter (mm)", "");
    const path = (ys) => xs.map((x, i) => `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(ys[i]).toFixed(1)}`).join("");
    // shade the top 2.3% (beyond z = +2) under the slider's shape
    const c2 = dAt(r, 2, state.lam);
    const tail = xs.map((x, i) => [x, dSel[i]]).filter(([x]) => x >= c2);
    if (tail.length > 1) svg.appendChild(el("path", { d: `M${X(c2)},${Y(0)}` + tail.map(([x, y]) => `L${X(x).toFixed(1)},${Y(y).toFixed(1)}`).join("") + `L${X(tail[tail.length - 1][0])},${Y(0)}Z`, style: `fill:var(${INFO[state.key].color});opacity:.2` }));
    svg.appendChild(el("path", { d: path(dFit), style: "fill:none;stroke:var(--muted);stroke-width:1.5;stroke-dasharray:5 4" }));
    svg.appendChild(el("path", { d: path(dSel), style: `fill:none;stroke:var(${INFO[state.key].color});stroke-width:2.5` }));
    const med = medMm(r);
    svg.appendChild(el("line", { x1: X(med), x2: X(med), y1: H - m.b, y2: H - m.b - 8, style: "stroke:var(--ink);stroke-width:2" }));
    svg.appendChild(el("text", { x: X(med), y: H - m.b - 12, "text-anchor": "middle", style: "fill:var(--ink);font:10.5px var(--font-num)" }, "median"));
    if (D <= x1) {
      svg.appendChild(el("line", { x1: X(D), x2: X(D), y1: m.t, y2: H - m.b, style: "stroke:var(--ink);stroke-width:1.5" }));
      svg.appendChild(el("text", { x: X(D) - 6, y: m.t - 6, "text-anchor": "end", style: "fill:var(--ink);font:600 11px var(--font-num)" }, D + " mm"));
    } else {
      svg.appendChild(el("text", { x: W - m.r, y: m.t - 6, "text-anchor": "end", style: "fill:var(--ink);font:600 11px var(--font-num)" }, `${D} mm is off this chart →`));
    }
    $("dens-legend").innerHTML = `<span><span class="dot" style="background:var(${INFO[state.key].color})"></span>λ ${state.lam.toFixed(2)} (slider)</span><span><span class="dot" style="background:var(--muted)"></span>λ ${fitted().toFixed(2)} (fitted, dashed)</span><span>Shaded: the top 2.3% (beyond z +2, from ${c2.toFixed(1)} mm) under the slider's shape</span>`;
  }

  function drawZ(r, D) {
    const svg = $("zplot"); while (svg.firstChild) svg.removeChild(svg.firstChild);
    const W = 640, H = 320, m = { l: 44, r: 112, t: 14, b: 44 };
    const med = medMm(r), x0 = Math.floor(med), x1 = Math.ceil(Math.max(D * 1.08, dAt(r, 6, 0)));
    const lams = LAMS();
    const zmax = Math.max(6, ...lams.map(([l]) => zOf(x1, r, l)));
    const y1 = Math.ceil(zmax / 2) * 2;
    const X = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), Y = (v) => m.t + (1 - v / y1) * (H - m.t - m.b);
    axes(svg, W, H, m, X, Y, ticks(x0, x1, x1 - x0 > 30 ? 10 : 5), ticks(0, y1, y1 > 12 ? 4 : 2), "Diameter (mm)", "z-score");
    [[2, "--band-mild"], [3.5, "--band-moderate"], [5, "--band-severe"]].forEach(([z, c]) => {
      svg.appendChild(el("line", { x1: m.l, x2: W - m.r, y1: Y(z), y2: Y(z), style: `stroke:var(${c});stroke-width:3;opacity:.9` }));
    });
    const xs = []; for (let i = 0; i <= 200; i++) xs.push(x0 + (x1 - x0) * i / 200);
    const dash = ["", "", "6 4", "2 4"];
    const ends = [];
    lams.forEach(([l, lab], j) => {
      const isFit = j === 1;
      svg.appendChild(el("path", { d: xs.map((x, i) => `${i ? "L" : "M"}${X(x).toFixed(1)},${Y(Math.min(zOf(x, r, l), y1)).toFixed(1)}`).join(""),
        style: `fill:none;stroke:${isFit ? `var(${INFO[state.key].color})` : "var(--ink)"};stroke-width:${isFit ? 2.5 : 1.5};opacity:${isFit ? 1 : .7};stroke-dasharray:${dash[j]}` }));
      ends.push({ y: Y(Math.min(zOf(x1, r, l), y1)), lab });
    });
    ends.sort((a, b) => a.y - b.y);
    for (let i = 1; i < ends.length; i++) if (ends[i].y - ends[i - 1].y < 13) ends[i].y = ends[i - 1].y + 13;
    ends.forEach((e) => svg.appendChild(el("text", { x: W - m.r + 6, y: e.y + 4, style: "fill:var(--ink);font:10.5px var(--font-num)" }, e.lab)));
    svg.appendChild(el("line", { x1: X(D), x2: X(D), y1: m.t, y2: H - m.b, style: "stroke:var(--ink);stroke-width:1.2;stroke-dasharray:3 3" }));
    lams.forEach(([l], j) => {
      const z = zOf(D, r, l);
      if (z <= y1) svg.appendChild(el("circle", { cx: X(D), cy: Y(z), r: 4.5, style: `fill:${j === 1 ? `var(${INFO[state.key].color})` : "var(--ink)"};stroke:var(--surface);stroke-width:2` }));
    });
    const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, style: "fill:transparent" });
    hit.addEventListener("mousemove", (e) => {
      const bb = svg.getBoundingClientRect(), sx = (e.clientX - bb.left) * W / bb.width;
      const d = x0 + (sx - m.l) / (W - m.l - m.r) * (x1 - x0);
      showTip(e, `<b>${d.toFixed(1)} mm</b> (${Math.round(100 * d / med)}% of expected)<br>` + lams.map(([l, lab]) => `${lab}: <b>${fmtZ(zOf(d, r, l))}</b>`).join("<br>"));
    });
    hit.addEventListener("mouseleave", () => { tip.hidden = true; });
    svg.appendChild(hit);
  }

  function renderTable(r, D) {
    const tb = $("lam-tab"); tb.textContent = "";
    const base = LAMS(), dup = base.some(([l]) => Math.abs(l - state.lam) < 0.005);
    (dup ? base : base.concat([[state.lam, `λ ${state.lam.toFixed(2)} (slider)`]])).forEach(([l, lab], j) => {
      const z = zOf(D, r, l), [cls, cat] = category(z);
      const tr = document.createElement("tr"); if (j === 1) tr.className = "primary";
      tr.innerHTML = `<td>${lab}</td><td class="num">${fmtZ(z)}</td><td><span class="chip ${cls}" style="font-size:12px">${cat}</span></td>` +
        `<td class="num">${dAt(r, 2, l).toFixed(1)}</td><td class="num">${dAt(r, 3.5, l).toFixed(1)}</td><td class="num">${dAt(r, 5, l).toFixed(1)}</td>`;
      tb.appendChild(tr);
    });
  }

  function update() {
    if (state.lam == null) state.lam = fitted();
    $("lam").value = state.lam; $("lam-out").textContent = state.lam.toFixed(2);
    $("fit-btn").textContent = `Fitted (${fitted().toFixed(2)})`;
    const c = child(), D = num("d");
    const ready = c.ok && D > 0;
    $("results").hidden = !ready; $("need").hidden = ready;
    if (!ready) return;
    const r = ref(state.key, c.p, c.hasWt);
    const z = zOf(D, r, state.lam), [cls, cat] = category(z);
    $("z").textContent = fmtZ(z); $("cat").className = "chip " + cls; $("cat").textContent = cat;
    $("pct").textContent = Math.round(100 * D / medMm(r)) + "%";
    $("med").textContent = medMm(r).toFixed(1) + " mm";
    $("sig").textContent = (100 * r.sigma).toFixed(1) + "%";
    $("model").textContent = c.hasWt ? "full model" : "height-and-sex model";
    drawDensity(r, D); drawZ(r, D); renderTable(r, D);
  }

  document.querySelectorAll("#level-tabs button").forEach((b) => b.addEventListener("click", () => {
    state.key = b.dataset.k; state.lam = null;
    document.querySelectorAll("#level-tabs button").forEach((x) => x.setAttribute("aria-pressed", x === b ? "true" : "false"));
    update();
  }));
  $("lam").addEventListener("input", () => { state.lam = +$("lam").value; update(); });
  document.querySelectorAll("[data-lam]").forEach((b) => b.addEventListener("click", () => {
    state.lam = b.dataset.lam === "fit" ? fitted() : +b.dataset.lam; update();
  }));
  $("f").addEventListener("input", (e) => { if (e.target.id !== "lam") update(); });
  $("f").addEventListener("submit", (e) => e.preventDefault());
  update();
})();
