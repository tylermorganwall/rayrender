(function () {
  const state = {
    rows: [],
    filtered: []
  };

  const metricConfigs = [
    { id: "renderChart", field: "render_seconds", label: "render_seconds" },
    { id: "bvhChart", field: "bvh_build_seconds", label: "bvh_build_seconds" },
    { id: "memoryChart", field: "max_rss_mb", label: "max_rss_mb" },
    {
      id: "packageChart",
      field: "package_build_seconds",
      fallback: "build_seconds",
      label: "package build seconds"
    }
  ];

  function embeddedCsv() {
    const node = document.getElementById("benchmark-csv");
    return node ? node.textContent.trim() : "";
  }

  function parseCsv(text) {
    const rows = [];
    let row = [];
    let field = "";
    let inQuotes = false;
    for (let i = 0; i < text.length; i += 1) {
      const ch = text[i];
      const next = text[i + 1];
      if (inQuotes) {
        if (ch === '"' && next === '"') {
          field += '"';
          i += 1;
        } else if (ch === '"') {
          inQuotes = false;
        } else {
          field += ch;
        }
      } else if (ch === '"') {
        inQuotes = true;
      } else if (ch === ",") {
        row.push(field);
        field = "";
      } else if (ch === "\n") {
        row.push(field);
        rows.push(row);
        row = [];
        field = "";
      } else if (ch !== "\r") {
        field += ch;
      }
    }
    if (field.length || row.length) {
      row.push(field);
      rows.push(row);
    }
    if (!rows.length) {
      return [];
    }
    const header = rows.shift();
    return rows
      .filter((item) => item.length && item.some((value) => value.length))
      .map((item) => {
        const out = {};
        header.forEach((name, index) => {
          out[name] = item[index] || "";
        });
        return out;
      });
  }

  function numberValue(row, field, fallback) {
    const value = row[field] || (fallback ? row[fallback] : "");
    if (!value || value === "NA") {
      return NaN;
    }
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : NaN;
  }

  function uniqueValues(rows, field) {
    return Array.from(new Set(rows.map((row) => row[field] || "NA"))).sort();
  }

  function populateSelect(id, values) {
    const select = document.getElementById(id);
    const previous = select.value;
    select.innerHTML = "";
    const all = document.createElement("option");
    all.value = "";
    all.textContent = "All";
    select.appendChild(all);
    values.forEach((value) => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      select.appendChild(option);
    });
    if (values.includes(previous)) {
      select.value = previous;
    }
  }

  function median(values) {
    const clean = values
      .map(Number)
      .filter((value) => Number.isFinite(value))
      .sort((a, b) => a - b);
    if (!clean.length) {
      return NaN;
    }
    const middle = Math.floor(clean.length / 2);
    return clean.length % 2 ? clean[middle] : (clean[middle - 1] + clean[middle]) / 2;
  }

  function groupKey(row) {
    return [
      row.branch_name || "NA",
      row.commit_sha || "NA",
      row.timestamp_utc || "NA",
      row.benchmark_name || "NA",
      row.build_config_name || "NA",
      row.width || "NA",
      row.height || "NA",
      row.samples || "NA",
      row.threads || "NA"
    ].join("\u001f");
  }

  function summarize(rows) {
    const groups = new Map();
    rows.forEach((row) => {
      const key = groupKey(row);
      if (!groups.has(key)) {
        groups.set(key, []);
      }
      groups.get(key).push(row);
    });
    return Array.from(groups.values()).map((items) => {
      const first = items[0];
      return {
        branch_name: first.branch_name || "NA",
        commit_sha: first.commit_sha || "NA",
        short_sha: (first.commit_sha || "NA").slice(0, 8),
        timestamp_utc: first.timestamp_utc || "NA",
        benchmark_name: first.benchmark_name || "NA",
        build_config_name: first.build_config_name || "NA",
        width: first.width || "NA",
        height: first.height || "NA",
        samples: first.samples || "NA",
        threads: first.threads || "NA",
        n: items.length,
        median_render_seconds: median(items.map((row) => numberValue(row, "render_seconds"))),
        min_render_seconds: Math.min(...items.map((row) => numberValue(row, "render_seconds")).filter(Number.isFinite)),
        max_render_seconds: Math.max(...items.map((row) => numberValue(row, "render_seconds")).filter(Number.isFinite)),
        median_bvh_build_seconds: median(items.map((row) => numberValue(row, "bvh_build_seconds"))),
        median_max_rss_mb: median(items.map((row) => numberValue(row, "max_rss_mb"))),
        median_total_seconds: median(items.map((row) => numberValue(row, "total_seconds"))),
        median_package_build_seconds: median(
          items.map((row) => numberValue(row, "package_build_seconds", "build_seconds"))
        )
      };
    });
  }

  function formatNumber(value) {
    return Number.isFinite(value) ? value.toFixed(3) : "NA";
  }

  function applyFilters() {
    const benchmark = document.getElementById("benchmarkFilter").value;
    const config = document.getElementById("configFilter").value;
    const branch = document.getElementById("branchFilter").value;
    state.filtered = state.rows.filter((row) => {
      return (!benchmark || row.benchmark_name === benchmark) &&
        (!config || row.build_config_name === config) &&
        (!branch || row.branch_name === branch);
    });
    render();
  }

  function latestRows(summary) {
    const sorted = [...summary].sort((a, b) => {
      return String(a.timestamp_utc).localeCompare(String(b.timestamp_utc));
    });
    const latest = new Map();
    sorted.forEach((row) => {
      latest.set([row.branch_name, row.benchmark_name, row.build_config_name].join("\u001f"), row);
    });
    return Array.from(latest.values()).sort((a, b) => {
      return [a.branch_name, a.benchmark_name, a.build_config_name]
        .join("\u001f")
        .localeCompare([b.branch_name, b.benchmark_name, b.build_config_name].join("\u001f"));
    });
  }

  function renderTable(summary) {
    const tbody = document.querySelector("#latestTable tbody");
    tbody.innerHTML = "";
    latestRows(summary).forEach((row) => {
      const tr = document.createElement("tr");
      [
        row.branch_name,
        row.benchmark_name,
        row.build_config_name,
        row.short_sha,
        row.timestamp_utc,
        row.n,
        formatNumber(row.median_render_seconds),
        formatNumber(row.median_bvh_build_seconds),
        formatNumber(row.median_max_rss_mb),
        formatNumber(row.median_total_seconds)
      ].forEach((value) => {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
  }

  function colorFor(value) {
    const palette = ["#2563eb", "#c2410c", "#047857", "#7c3aed", "#be123c", "#0f766e"];
    let hash = 0;
    for (let i = 0; i < value.length; i += 1) {
      hash = (hash * 31 + value.charCodeAt(i)) >>> 0;
    }
    return palette[hash % palette.length];
  }

  function renderChart(config, summary, rows) {
    const host = document.getElementById(config.id);
    host.innerHTML = "";
    const values = summary
      .map((row) => ({
        row,
        x: Date.parse(row.timestamp_utc),
        y: config.field === "render_seconds" ? row.median_render_seconds :
          config.field === "bvh_build_seconds" ? row.median_bvh_build_seconds :
          config.field === "max_rss_mb" ? row.median_max_rss_mb :
          row.median_package_build_seconds
      }))
      .filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y));
    if (!values.length) {
      host.textContent = "not available";
      host.className = "chart empty";
      return;
    }
    host.className = "chart";
    const width = 900;
    const height = 260;
    const margin = { top: 18, right: 18, bottom: 44, left: 58 };
    const minX = Math.min(...values.map((point) => point.x));
    const maxX = Math.max(...values.map((point) => point.x));
    const minY = Math.min(0, Math.min(...values.map((point) => point.y)));
    const maxY = Math.max(...values.map((point) => point.y));
    const spanX = Math.max(1, maxX - minX);
    const spanY = Math.max(1e-9, maxY - minY);
    const sx = (value) => margin.left + ((value - minX) / spanX) * (width - margin.left - margin.right);
    const sy = (value) => height - margin.bottom - ((value - minY) / spanY) * (height - margin.top - margin.bottom);
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", config.label);

    const axis = document.createElementNS(svg.namespaceURI, "path");
    axis.setAttribute("d", `M${margin.left},${margin.top}V${height - margin.bottom}H${width - margin.right}`);
    axis.setAttribute("class", "axis");
    svg.appendChild(axis);

    const yLabel = document.createElementNS(svg.namespaceURI, "text");
    yLabel.setAttribute("x", 8);
    yLabel.setAttribute("y", 18);
    yLabel.setAttribute("class", "axis-label");
    yLabel.textContent = config.label;
    svg.appendChild(yLabel);

    const grouped = new Map();
    values.forEach((point) => {
      const key = `${point.row.branch_name} / ${point.row.build_config_name}`;
      if (!grouped.has(key)) {
        grouped.set(key, []);
      }
      grouped.get(key).push(point);
    });
    grouped.forEach((points, key) => {
      points.sort((a, b) => a.x - b.x);
      const path = document.createElementNS(svg.namespaceURI, "path");
      path.setAttribute(
        "d",
        points.map((point, index) => `${index ? "L" : "M"}${sx(point.x)},${sy(point.y)}`).join("")
      );
      path.setAttribute("fill", "none");
      path.setAttribute("stroke", colorFor(key));
      path.setAttribute("stroke-width", "2");
      svg.appendChild(path);
    });

    rows.forEach((row) => {
      const x = Date.parse(row.timestamp_utc);
      const y = numberValue(row, config.field, config.fallback);
      if (!Number.isFinite(x) || !Number.isFinite(y)) {
        return;
      }
      const point = document.createElementNS(svg.namespaceURI, "circle");
      point.setAttribute("cx", sx(x));
      point.setAttribute("cy", sy(y));
      point.setAttribute("r", "2.5");
      point.setAttribute("class", "raw-point");
      svg.appendChild(point);
    });

    values.forEach((point) => {
      const circle = document.createElementNS(svg.namespaceURI, "circle");
      circle.setAttribute("cx", sx(point.x));
      circle.setAttribute("cy", sy(point.y));
      circle.setAttribute("r", "4");
      circle.setAttribute("fill", colorFor(`${point.row.branch_name} / ${point.row.build_config_name}`));
      const title = document.createElementNS(svg.namespaceURI, "title");
      title.textContent = `${point.row.branch_name} ${point.row.build_config_name} ${point.row.short_sha}: ${formatNumber(point.y)}`;
      circle.appendChild(title);
      svg.appendChild(circle);
    });
    host.appendChild(svg);
  }

  function render() {
    const summary = summarize(state.filtered);
    document.getElementById("rowCount").textContent = `${state.filtered.length} rows`;
    renderTable(summary);
    metricConfigs.forEach((config) => renderChart(config, summary, state.filtered));
  }

  async function loadData() {
    let text = embeddedCsv();
    if (!text) {
      const response = await fetch("data/render_benchmarks.csv");
      text = await response.text();
    }
    state.rows = parseCsv(text);
    state.filtered = state.rows;
    populateSelect("benchmarkFilter", uniqueValues(state.rows, "benchmark_name"));
    populateSelect("configFilter", uniqueValues(state.rows, "build_config_name"));
    populateSelect("branchFilter", uniqueValues(state.rows, "branch_name"));
    ["benchmarkFilter", "configFilter", "branchFilter"].forEach((id) => {
      document.getElementById(id).addEventListener("change", applyFilters);
    });
    render();
  }

  loadData().catch((error) => {
    document.getElementById("rowCount").textContent = `failed to load data: ${error.message}`;
  });
}());
