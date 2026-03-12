#!/usr/bin/env Rscript

# BaitBench Cross-Reactivity Report Generator
# Renders an HTML report from xreact analysis results

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--hits", type = "character", help = "Hits TSV file"),
  make_option("--summary", type = "character", help = "Summary TSV file"),
  make_option("--params", type = "character", default = NULL, help = "Run parameters TSV file"),
  make_option("--threshold", type = "double", default = 80.0, help = "Homology threshold [default %default]"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$hits) || is.null(opt$summary) || is.null(opt$output)) {
  stop("--hits, --summary, and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "xreact.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find xreact.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    hits_file = normalizePath(opt$hits),
    summary_file = normalizePath(opt$summary),
    params_file = if (!is.null(opt$params)) normalizePath(opt$params) else "",
    threshold = opt$threshold
  ),
  quiet = TRUE
)

message("Cross-reactivity report generated: ", opt$output)
