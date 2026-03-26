#!/usr/bin/env Rscript

# BaitBench Build Probes Report Generator
# Renders an HTML report from probe building pipeline statistics

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--stats", type = "character", help = "Pipeline statistics TSV file"),
  make_option("--params", type = "character", default = NULL, help = "Run parameters TSV file"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$stats) || is.null(opt$output)) {
  stop("--stats and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "build_probes.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find build_probes.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    stats_file = normalizePath(opt$stats),
    params_file = if (!is.null(opt$params)) normalizePath(opt$params) else ""
  ),
  quiet = TRUE
)

message("Build probes report generated: ", opt$output)
