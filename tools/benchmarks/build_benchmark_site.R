script = sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
source(file.path(dirname(script), "benchmark_data.R"))
args = benchmark_cli(c("history_csv", "site_dir"))
# Missing/corrupt established history is a deployment error, never an empty site.
rows = benchmark_read(args$history_csv)
data = benchmark_dashboard(rows)
assets = file.path(dirname(script), "site_assets")
dir.create(
  file.path(args$site_dir, "assets"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(args$site_dir, "data"),
  recursive = TRUE,
  showWarnings = FALSE
)
stopifnot(file.copy(
  args$history_csv,
  file.path(args$site_dir, "data", "render_benchmarks.csv"),
  overwrite = TRUE
))
for (name in c("dashboard.js", "style.css")) {
  stopifnot(file.copy(
    file.path(assets, name),
    file.path(args$site_dir, "assets", name),
    overwrite = TRUE
  ))
}
json = benchmark_json(data)
writeLines(json, file.path(args$site_dir, "data", "dashboard.json"))
# A script element is raw text: HTML entities are not decoded. JSON unicode
# escapes preserve exact strings and prevent history text closing the element.
json = gsub("<", "\\u003c", json, fixed = TRUE)
html = readLines(file.path(assets, "index.html"), warn = FALSE)
html[html == "BENCHMARK_DATA"] = json
writeLines(html, file.path(args$site_dir, "index.html"))
invisible(file.create(file.path(args$site_dir, ".nojekyll")))
message("Benchmark site written to ", args$site_dir)
