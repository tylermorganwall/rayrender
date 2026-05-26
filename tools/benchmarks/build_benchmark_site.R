usage = function() {
  cat(
    paste(
      "Usage:",
      "Rscript tools/benchmarks/build_benchmark_site.R \\",
      "  --history-csv data/render_benchmarks.csv \\",
      "  --site-dir site",
      sep = "\n"
    )
  )
}

stop_usage = function(message) {
  usage()
  stop(message, call. = FALSE)
}

parse_args = function(argv) {
  values = list()
  i = 1
  while (i <= length(argv)) {
    key = argv[[i]]
    if (key %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key, "--") || i == length(argv)) {
      stop_usage(paste0("Invalid argument: ", key))
    }
    name = gsub("-", "_", sub("^--", "", key))
    values[[name]] = argv[[i + 1]]
    i = i + 2
  }
  missing = setdiff(c("history_csv", "site_dir"), names(values))
  if (length(missing) > 0) {
    stop_usage(paste0("Missing arguments: ", paste(missing, collapse = ", ")))
  }
  values
}

script_dir = function() {
  file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  getwd()
}

html_escape = function(value) {
  value = gsub("&", "&amp;", value, fixed = TRUE)
  value = gsub("<", "&lt;", value, fixed = TRUE)
  value = gsub(">", "&gt;", value, fixed = TRUE)
  value
}

copy_asset = function(from, to) {
  if (!file.exists(from)) {
    stop("Missing dashboard asset: ", from, call. = FALSE)
  }
  invisible(file.copy(from, to, overwrite = TRUE))
}

write_latest_json = function(path) {
  if (file.exists(path)) {
    return(invisible(TRUE))
  }
  writeLines('{"row_count":0,"latest":[]}', path)
}

args = parse_args(commandArgs(trailingOnly = TRUE))
history_csv = normalizePath(args$history_csv, mustWork = TRUE)
site_dir = normalizePath(args$site_dir, mustWork = FALSE)
site_data_dir = file.path(site_dir, "data")
site_asset_dir = file.path(site_dir, "assets")
source_asset_dir = file.path(script_dir(), "site_assets")

dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(site_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(site_asset_dir, recursive = TRUE, showWarnings = FALSE)

site_csv = file.path(site_data_dir, "render_benchmarks.csv")
site_latest = file.path(site_data_dir, "latest.json")

invisible(file.copy(history_csv, site_csv, overwrite = TRUE))
latest_source = file.path(dirname(history_csv), "latest.json")
if (file.exists(latest_source)) {
  invisible(file.copy(latest_source, site_latest, overwrite = TRUE))
} else {
  write_latest_json(site_latest)
}

copy_asset(
  file.path(source_asset_dir, "dashboard.js"),
  file.path(site_asset_dir, "dashboard.js")
)
copy_asset(
  file.path(source_asset_dir, "style.css"),
  file.path(site_asset_dir, "style.css")
)
invisible(file.create(file.path(site_dir, ".nojekyll")))

csv_text = paste(readLines(history_csv, warn = FALSE), collapse = "\n")
index = paste(
  '<!doctype html>',
  '<html lang="en">',
  '<head>',
  '  <meta charset="utf-8">',
  '  <meta name="viewport" content="width=device-width, initial-scale=1">',
  '  <title>rayrender Benchmarks</title>',
  '  <link rel="stylesheet" href="assets/style.css">',
  '</head>',
  '<body>',
  '  <header>',
  '    <h1>rayrender render benchmarks</h1>',
  '    <p>Persistent benchmark history from the benchmark-history branch.</p>',
  '  </header>',
  '  <main>',
  '    <section class="toolbar" aria-label="Filters">',
  '      <div class="field">',
  '        <label for="benchmarkFilter">Benchmark</label>',
  '        <select id="benchmarkFilter"></select>',
  '      </div>',
  '      <div class="field">',
  '        <label for="configFilter">Build config</label>',
  '        <select id="configFilter"></select>',
  '      </div>',
  '      <div class="field">',
  '        <label for="branchFilter">Branch</label>',
  '        <select id="branchFilter"></select>',
  '      </div>',
  '      <div class="meta">',
  '        <span id="rowCount">loading</span> &middot; ',
  '        <a href="data/render_benchmarks.csv">raw CSV</a>',
  '      </div>',
  '    </section>',
  '    <section>',
  '      <h2>Latest medians</h2>',
  '      <div class="table-wrap">',
  '        <table id="latestTable">',
  '          <thead>',
  '            <tr>',
  '              <th>Branch</th>',
  '              <th>Benchmark</th>',
  '              <th>Config</th>',
  '              <th>Commit</th>',
  '              <th>Timestamp</th>',
  '              <th>n</th>',
  '              <th>Render median</th>',
  '              <th>Render mean</th>',
  '              <th>BVH build median</th>',
  '              <th>Max RSS median</th>',
  '              <th>Total median</th>',
  '            </tr>',
  '          </thead>',
  '          <tbody></tbody>',
  '        </table>',
  '      </div>',
  '    </section>',
  '    <section class="charts">',
  '      <div class="chart-panel"><h2>Render seconds over time</h2><div id="renderChart" class="chart"></div></div>',
  '      <div class="chart-panel"><h2>BVH build seconds over time</h2><div id="bvhChart" class="chart"></div></div>',
  '      <div class="chart-panel"><h2>Max RSS MB over time</h2><div id="memoryChart" class="chart"></div></div>',
  '      <div class="chart-panel"><h2>Package build seconds over time</h2><div id="packageChart" class="chart"></div></div>',
  '    </section>',
  '  </main>',
  '  <script type="text/plain" id="benchmark-csv">',
  html_escape(csv_text),
  '  </script>',
  '  <script src="assets/dashboard.js"></script>',
  '</body>',
  '</html>',
  sep = "\n"
)

writeLines(index, file.path(site_dir, "index.html"))
message("Benchmark site written to ", site_dir)
