#!/usr/bin/env Rscript

# BaitBench CT Sweep Report Generator
# Renders an HTML report from CT sweep depth curve data

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--sweep", type = "character", help = "CT sweep depth curves TSV"),
  make_option("--sample-ids", type = "character", help = "Comma-separated sample target IDs"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sweep) || is.null(opt$`sample-ids`) || is.null(opt$output)) {
  stop("--sweep, --sample-ids, and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "ct_sweep.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find ct_sweep.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    sweep_file = normalizePath(opt$sweep),
    sample_ids = strsplit(opt$`sample-ids`, ",")[[1]]
  ),
  quiet = TRUE
)

message("CT sweep report generated: ", opt$output)
