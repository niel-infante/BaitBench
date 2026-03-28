#!/usr/bin/env Rscript

# BaitBench Probe Assessment Report Generator
# Combined report: build pipeline (optional) + coverage + cross-reactivity

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--build-stats", type = "character", default = NULL,
              help = "Build pipeline statistics TSV (optional, from build-probes)"),
  make_option("--build-params", type = "character", default = NULL,
              help = "Build pipeline parameters TSV (optional, from build-probes)"),
  make_option("--xreact-hits", type = "character",
              help = "Cross-reactivity hits TSV"),
  make_option("--xreact-summary", type = "character",
              help = "Cross-reactivity summary TSV"),
  make_option("--threshold", type = "double", default = 80.0,
              help = "Homology threshold [default %default]"),
  make_option("--cov-summary", type = "character",
              help = "Probe coverage summary TSV"),
  make_option("--cov-depth", type = "character",
              help = "Probe depth intervals TSV"),
  make_option("--cov-multi-mapping", type = "character", default = NULL,
              help = "Multi-mapping probes TSV"),
  make_option("--proximity", type = "integer", default = 50,
              help = "Pull-down zone distance [default %default]"),
  make_option("--params", type = "character", default = NULL,
              help = "Overall run parameters TSV"),
  make_option("--output", type = "character",
              help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$`xreact-hits`) || is.null(opt$`xreact-summary`) ||
    is.null(opt$`cov-summary`) || is.null(opt$`cov-depth`) ||
    is.null(opt$output)) {
  stop("Required: --xreact-hits, --xreact-summary, --cov-summary, --cov-depth, --output")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "assess_probes.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find assess_probes.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    build_stats_file = if (!is.null(opt$`build-stats`)) normalizePath(opt$`build-stats`) else "",
    build_params_file = if (!is.null(opt$`build-params`)) normalizePath(opt$`build-params`) else "",
    xreact_hits_file = normalizePath(opt$`xreact-hits`),
    xreact_summary_file = normalizePath(opt$`xreact-summary`),
    xreact_threshold = opt$threshold,
    coverage_summary_file = normalizePath(opt$`cov-summary`),
    coverage_depth_file = normalizePath(opt$`cov-depth`),
    coverage_multi_mapping_file = if (!is.null(opt$`cov-multi-mapping`))
      normalizePath(opt$`cov-multi-mapping`) else "",
    coverage_proximity = opt$proximity,
    params_file = if (!is.null(opt$params)) normalizePath(opt$params) else ""
  ),
  quiet = TRUE
)

message("Probe assessment report generated: ", opt$output)
