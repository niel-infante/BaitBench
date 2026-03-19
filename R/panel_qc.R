#!/usr/bin/env Rscript

# BaitBench Panel QC Report Generator
# Renders an HTML report from target panel discriminability analysis

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--discriminability", type = "character", help = "Species discriminability TSV file"),
  make_option("--matrix", type = "character", help = "Species confusion matrix TSV file"),
  make_option("--similarity", type = "character", help = "Target similarity TSV file"),
  make_option("--params", type = "character", default = NULL, help = "Run parameters TSV file"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$discriminability) || is.null(opt$matrix) || is.null(opt$output)) {
  stop("--discriminability, --matrix, and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "panel_qc.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find panel_qc.Rmd in ", script_dir)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    discriminability_file = normalizePath(opt$discriminability),
    matrix_file = normalizePath(opt$matrix),
    similarity_file = if (!is.null(opt$similarity)) normalizePath(opt$similarity) else "",
    params_file = if (!is.null(opt$params)) normalizePath(opt$params) else ""
  ),
  quiet = TRUE
)

message("Panel QC report generated: ", opt$output)
