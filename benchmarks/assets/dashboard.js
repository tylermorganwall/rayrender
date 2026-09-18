(async function () {
  "use strict";
  const node = (id) => document.getElementById(id);
  const text = (id, value) => { node(id).textContent = value; };
  const available = (value) => typeof value === "number" && Number.isFinite(value);
  const numberText = (value) => value > 0 && value < 0.001 ? value.toPrecision(4) : value.toFixed(4);
  const label = (attempt) => `run ${attempt.run_id || "legacy"} / attempt ${attempt.run_attempt || "unknown"} / ${(attempt.commit_sha || "unknown").slice(0, 12)} · ${attempt.status}`;
  const metric = (row, name, unit) => available(row[name]) ? `${numberText(row[name])} ${unit} (n=${row[`${name}_n`]})` : name === "bvh_build_seconds" ? "BVH timing not recorded" : "Unavailable";
  const settings = (row) => `${row.benchmark_settings_json || `${row.width}×${row.height}, ${row.samples} samples, ${row.threads} threads, seed ${row.seed}`} · warmups ${row.warmup_iterations || "unrecorded"} · requested iterations ${row.iterations || "unrecorded"}`;
  const list = (value) => Array.isArray(value) ? value : value ? [value] : [];
  function select(id, values, all = true) {
    node(id).replaceChildren();
    if (all) node(id).add(new Option("All", ""));
    values.forEach(([value, title]) => node(id).add(new Option(title, value)));
  }
  function addText(parent, tag, value) {
    const child = document.createElement(tag);
    child.textContent = value;
    parent.appendChild(child);
    return child;
  }
  function chart(id, field, rows) {
    const host = node(id);
    host.replaceChildren();
    const points = rows.filter((r) => available(r[field]) && Number.isFinite(Date.parse(r.timestamp_utc)));
    if (!points.length) {
      host.className = "chart empty";
      host.textContent = field === "bvh_build_seconds" ? "BVH timing not recorded" : "No successful measurements available";
      return;
    }
    host.className = "chart";
    const ns = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(ns, "svg");
    svg.setAttribute("viewBox", "0 0 900 230");
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", field);
    const times = points.map((r) => Date.parse(r.timestamp_utc));
    const low = Math.min(...times), high = Math.max(...times);
    const max = Math.max(...points.map((r) => r[field]), 0.000001);
    function svgText(x, y, content) {
      const n = document.createElementNS(ns, "text"); n.setAttribute("x", x); n.setAttribute("y", y); n.textContent = content; svg.appendChild(n);
    }
    svgText(4, 20, max.toFixed(3)); svgText(4, 195, "0");
    svgText(60, 220, new Date(low).toISOString().slice(0, 10));
    svgText(770, 220, new Date(high).toISOString().slice(0, 10));
    points.forEach((r) => {
      const circle = document.createElementNS(ns, "circle");
      circle.setAttribute("cx", 70 + (Date.parse(r.timestamp_utc) - low) / Math.max(1, high - low) * 760);
      circle.setAttribute("cy", 195 - r[field] / max * 175);
      circle.setAttribute("r", 5);
      circle.setAttribute("fill", r.effective_backend === "scalar" ? "#c2410c" : "#2563eb");
      const title = document.createElementNS(ns, "title");
      title.textContent = `${r.build_config_name} revision ${r.config_revision || "legacy"}, ${r.benchmark_name}, run ${r.run_id}/${r.run_attempt}, ${r[field]} (n=${r[`${field}_n`]}), ${settings(r)}, ${r.effective_flags || "flags unrecorded"}`;
      circle.appendChild(title); svg.appendChild(circle);
    });
    host.appendChild(svg);
  }
  try {
    const embedded = node("benchmark-data").textContent.trim();
    let data;
    if (embedded) data = JSON.parse(embedded);
    else {
      const response = await fetch("data/dashboard.json");
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      data = await response.json();
    }
    if (!data || !Array.isArray(data.summaries) || !Array.isArray(data.attempts) || !Array.isArray(data.builds)) throw new Error("Invalid dashboard data schema");
    data.attempts.sort((a, b) => b.timestamp_utc.localeCompare(a.timestamp_utc));
    const values = (field) => [...new Set(data.summaries.map((r) => r[field] || "unknown"))].sort().map((x) => [x, x]);
    select("branchFilter", values("branch_name"));
    select("benchmarkFilter", values("benchmark_name"));
    select("configFilter", values("build_config_name"));
    // Start on the newest attempt's branch, and prefer its latest validated comparison.
    if (data.attempts.length) node("branchFilter").value = data.attempts[0].branch_name;
    function attempts() { return data.attempts.filter((a) => !node("branchFilter").value || a.branch_name === node("branchFilter").value); }
    function updateRuns() {
      const items = attempts();
      select("runFilter", items.map((a) => [a.attempt_key, label(a)]), false);
      const complete = items.find((a) => a.status === "complete");
      if (complete) node("runFilter").value = complete.attempt_key;
      text("latestAttempt", items.length ? `Latest attempt: ${label(items[0])} · ${items[0].timestamp_utc}` : "No benchmark attempts recorded");
      text("latestComplete", complete ? `Latest complete valid comparison: ${label(complete)}` : "No complete valid comparison recorded. Historical measurements remain available.");
    }
    function render() {
      const history = data.summaries.filter((r) => ["branch", "benchmark", "config"].every((key) => {
        const value = node(`${key}Filter`).value;
        const field = {branch: "branch_name", benchmark: "benchmark_name", config: "build_config_name"}[key];
        return !value || r[field] === value;
      }));
      const selected = history.filter((r) => r.attempt_key === node("runFilter").value);
      const attempt = data.attempts.find((a) => a.attempt_key === node("runFilter").value);
      text("rowCount", `${selected.reduce((n, r) => n + r.row_count, 0)} selected rows / ${data.row_count} historical rows`);
      text("selectionStatus", attempt ? `Displaying ${label(attempt)}. ${selected.length ? "" : "No measurements match these filters."}` : "No benchmark data available.");
      text("comparisonNote", new Set(selected.map(settings)).size > 1 ? "Different benchmark settings are shown separately; these rows are not equivalent comparisons." : "Configuration revisions and settings identify the measurements; legacy flags are unverified.");
      const body = document.querySelector("#latestTable tbody"); body.replaceChildren();
      selected.forEach((r) => {
        const tr = document.createElement("tr");
        const status = Object.entries(r.statuses).map(([k, n]) => `${k}: ${n}`).join(", ");
        const cells = [`${r.run_id || "legacy"} / ${r.run_attempt || "?"} / ${(r.commit_sha || "unknown").slice(0, 12)}`, r.benchmark_name, `${r.build_config_name} / rev ${r.config_revision || "legacy"} / ${r.effective_backend || "unrecorded"}`, `${r.width}×${r.height} · ${r.samples} samples · ${r.threads} threads · seed ${r.seed} · ${r.warmup_iterations || "unrecorded"} warmups`, status, metric(r, "render_seconds", "s"), metric(r, "bvh_build_seconds", "s"), metric(r, "max_rss_mb", "MB"), metric(r, "total_seconds", "s")];
        cells.forEach((v) => addText(tr, "td", v));
        const settingDetails = addText(tr.children[3], "details", "");
        addText(settingDetails, "summary", "All settings");
        addText(settingDetails, "pre", settings(r));
        const details = addText(tr.children[2], "details", "");
        addText(details, "summary", "Compiler and flags");
        addText(details, "p", `${r.cxx_standard || "Standard unrecorded"}; ${r.compiler_id || "Compiler unrecorded"}; ${r.effective_flags || "Effective flags unrecorded"}`);
        body.appendChild(tr);
      });
      const builds = node("buildMeasurements"); builds.replaceChildren();
      data.builds.filter((b) => b.attempt_key === node("runFilter").value && (!node("configFilter").value || b.build_config_name === node("configFilter").value)).forEach((b) => addText(builds, "p", `${b.build_config_name} / rev ${b.config_revision || "legacy"}: ${available(b.seconds) ? `${b.seconds.toFixed(3)} s (n=${b.n} build)` : "No successful build measurement"}`));
      const diagnostics = node("diagnostics"); diagnostics.replaceChildren();
      // Latest failures stay visible even while the table shows the last good comparison.
      const diagnosticKeys = new Set([node("runFilter").value, attempts()[0]?.attempt_key]);
      data.attempts.filter((a) => diagnosticKeys.has(a.attempt_key)).forEach((a) => {
        addText(diagnostics, "h3", label(a));
        list(a.validation_errors).forEach((error) => addText(diagnostics, "p", error));
        if (/^[0-9]+$/.test(a.run_id)) {
          const link = addText(diagnostics, "a", "Actions run and diagnostic artifacts");
          link.href = `https://github.com/tylermorganwall/rayrender/actions/runs/${a.run_id}/attempts/${a.run_attempt}`;
        }
        data.summaries.filter((r) => r.attempt_key === a.attempt_key).forEach((r) => {
          const failures = list(r.failures);
          if (failures.length) {
            const details = addText(diagnostics, "details", "");
            addText(details, "summary", `${r.build_config_name} / ${r.benchmark_name}: ${Object.keys(r.statuses).join(", ")}`);
            failures.forEach((error) => addText(details, "pre", error));
          }
        });
      });
      chart("renderChart", "render_seconds", history); chart("bvhChart", "bvh_build_seconds", history); chart("memoryChart", "max_rss_mb", history);
    }
    updateRuns(); render();
    node("branchFilter").addEventListener("change", () => { updateRuns(); render(); });
    ["runFilter", "benchmarkFilter", "configFilter"].forEach((id) => node(id).addEventListener("change", render));
    text("loadStatus", data.row_count ? "Benchmark data loaded." : "No benchmark data available.");
    node("loadStatus").dataset.state = data.row_count ? "loaded" : "empty";
  } catch (error) {
    text("loadStatus", `Dashboard load error: ${error.message}`);
    node("loadStatus").dataset.state = "error";
  }
}());
