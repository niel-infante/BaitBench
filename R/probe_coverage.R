#!/usr/bin/env Rscript

# BaitBench Probe Coverage Report Generator
# Renders an HTML report from probe coverage analysis results

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--summary", type = "character", help = "Probe coverage summary TSV"),
  make_option("--depth", type = "character", help = "Per-position probe depth TSV"),
  make_option("--multi-mapping", type = "character", default = NULL, help = "Multi-mapping probes TSV"),
  make_option("--proximity", type = "integer", default = 50, help = "Pull-down zone distance in bp [default %default]"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$summary) || is.null(opt$depth) || is.null(opt$output)) {
  stop("--summary, --depth, and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "probe_coverage.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find probe_coverage.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    summary_file = normalizePath(opt$summary),
    depth_file = normalizePath(opt$depth),
    multi_mapping_file = if (!is.null(opt$`multi-mapping`)) normalizePath(opt$`multi-mapping`) else "",
    proximity = opt$proximity
  ),
  quiet = TRUE
)

message("Probe coverage report generated: ", opt$output)
