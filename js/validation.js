/* Validation page: recomputes every reference case in the browser and compares with R. */
(function () {
  "use strict";
  const { INFO, selfTest, M } = AorticModel;
  const $ = (id) => document.getElementById(id);
  const res = selfTest();
  const z3 = (z) => (z >= 0 ? "+" : "−") + Math.abs(z).toFixed(3);
  const fmtD = (d) => (d >= 0 ? "+" : "−") + Math.abs(d).toExponential(1);

  $("s-n").textContent = res.n;
  $("s-z").textContent = 2 * res.n;
  $("s-max").textContent = res.maxd.toExponential(1);
  const st = $("s-status");
  st.textContent = res.pass ? "All match" : "Mismatch";
  st.className = "pill " + (res.pass ? "ok" : "bad");
  $("s-date").textContent = M.meta.created;

  let filter = "ALL";
  function render() {
    const tb = $("cases"); tb.textContent = "";
    res.rows.filter((r) => filter === "ALL" || r.structure === filter).forEach((r) => {
      const ok = Math.abs(r.diff) < 0.001 && Math.abs(r.diff_hs) < 0.001;
      const tr = document.createElement("tr");
      tr.innerHTML =
        `<td style="white-space:nowrap"><span class="dot" style="background:var(${INFO[r.structure].color})"></span>${INFO[r.structure].short}</td>` +
        `<td>${r.sex === "M" ? "Male" : "Female"}</td>` +
        `<td class="num">${r.age_y.toFixed(2)}</td><td class="num">${r.height_cm.toFixed(1)}</td>` +
        `<td class="num">${r.weight_kg.toFixed(1)}</td><td class="num">${r.diam_mm}</td>` +
        `<td class="num">${z3(r.z_R)}</td><td class="num">${z3(r.z_js)}</td><td class="num sub">${fmtD(r.diff)}</td>` +
        `<td class="num">${z3(r.z_R_hs)}</td><td class="num">${z3(r.z_js_hs)}</td><td class="num sub">${fmtD(r.diff_hs)}</td>` +
        `<td><span class="pill ${ok ? "ok" : "bad"}">${ok ? "Match" : "Check"}</span></td>`;
      tb.appendChild(tr);
    });
  }
  document.querySelectorAll("#filter button").forEach((btn) => btn.addEventListener("click", () => {
    filter = btn.dataset.k;
    document.querySelectorAll("#filter button").forEach((b) => b.setAttribute("aria-pressed", b === btn ? "true" : "false"));
    render();
  }));
  render();

  // ---------- Lauzon-Schnittka check ----------
  if (window.LauzonSchnittka && window.LS_TESTS) {
    const ls = LauzonSchnittka.selfTest();
    $("ls-sec").hidden = false;
    $("ls-sum").innerHTML = `<span class="pill ${ls.pass ? "ok" : "bad"}">${ls.pass ? "All match" : "Mismatch"}</span> ${3 * ls.n} z-scores and ${ls.n} BMI z-scores; largest difference ${Math.max(ls.maxd, ls.maxb).toExponential(1)}.`;
    const z3 = (z) => (z >= 0 ? "+" : "−") + Math.abs(z).toFixed(3);
    ls.rows.forEach(({ t, z, bmiz }) => {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td>${t.sex === "m" ? "Male" : "Female"}</td><td class="num">${t.age_y.toFixed(2)}</td><td class="num">${t.ht}</td><td class="num">${t.wt}</td>` +
        `<td class="num">${z3(t.bmiz)}</td><td class="num">${z3(bmiz)}</td>` +
        ["SOV", "STJ", "AAO"].map((k) => `<td class="num">${z3(t.z[k])}</td><td class="num">${z3(z[k])}</td>`).join("");
      $("ls-cases").appendChild(tr);
    });
  }
})();
