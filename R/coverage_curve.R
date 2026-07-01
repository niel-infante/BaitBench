#!/usr/bin/env Rscript

# BaitBench Coverage Curve Report Generator
# Renders an HTML report from coverage curve depth data

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--sweep", type = "character", help = "Coverage curve depth curves TSV"),
  make_option("--summary", type = "character", default = NULL, help = "Coverage curve summary stats TSV"),
  make_option("--sample-ids", type = "character", help = "Comma-separated sample target IDs"),
  make_option("--swept-params", type = "character", default = "",
              help = "Comma-separated swept parameter names (ct, fold_enrichment, num_sequences)"),
  make_option("--params", type = "character", default = NULL, help = "Run parameters TSV file"),
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

rmd_path <- file.path(script_dir, "coverage_curve.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find coverage_curve.Rmd in ", script_dir)
}

swept_params <- if (nchar(opt$`swept-params`) > 0) {
  strsplit(opt$`swept-params`, ",")[[1]]
} else {
  character(0)
}

rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    sweep_file = normalizePath(opt$sweep),
    summary_file = if (!is.null(opt$summary)) normalizePath(opt$summary) else "",
    sample_ids = strsplit(opt$`sample-ids`, ",")[[1]],
    swept_params = swept_params,
    params_file = if (!is.null(opt$params)) normalizePath(opt$params) else ""
  ),
  quiet = TRUE
)

message("Coverage curve report generated: ", opt$output)
