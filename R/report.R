#!/usr/bin/env Rscript

# BaitBench Report Generator
# Renders an HTML report from pipeline results using RMarkdown

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--summary", type = "character", help = "Summary TSV file (results.tsv)"),
  make_option("--detail", type = "character", help = "Detail TSV file (detected_detail.tsv)"),
  make_option("--run-name", type = "character", default = "BaitBench Run", help = "Run name"),
  make_option("--output", type = "character", help = "Output HTML file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$summary) || is.null(opt$detail) || is.null(opt$output)) {
  stop("--summary, --detail, and --output are required")
}

# Find the Rmd template relative to this script
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("--file=", commandArgs(trailingOnly = FALSE))
] |> sub("--file=", "", x = _)))

rmd_path <- file.path(script_dir, "report.Rmd")
if (!file.exists(rmd_path)) {
  stop("Cannot find report.Rmd in ", script_dir)
}

# Render the report
rmarkdown::render(
  input = rmd_path,
  output_file = normalizePath(opt$output, mustWork = FALSE),
  params = list(
    summary_file = normalizePath(opt$summary),
    detail_file = normalizePath(opt$detail),
    run_name = opt$`run-name`
  ),
  quiet = TRUE
)

message("Report generated: ", opt$output)
