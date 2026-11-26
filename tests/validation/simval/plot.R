#!/usr/bin/env Rscript
# Scatter of real vs simulated per-organism composition for the D6331 mock community.
#
#   usage: plot.R <comparison.tsv> <output.png> [title]
#
# The TSV comes from `SIMVAL_TSV=... compare.py <tag>`. Proportions span ~5 orders of
# magnitude, so both axes are log10; perfect agreement is the diagonal, and the shaded
# band is 2-fold either side. Organisms absent from one side are pinned to a floor and
# drawn hollow, because a zero cannot be placed on a log axis and should not be mistaken
# for a measurement.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot.R <comparison.tsv> <output.png> [title]")
tsv <- args[1]; out <- args[2]
title <- if (length(args) >= 3) args[3] else "D6331 mock community: simulated vs real composition"

d <- read.delim(tsv, stringsAsFactors = FALSE)
d <- d[d$set == "E. coli strains collapsed", ]
if (!nrow(d)) stop("no rows for the collapsed set in ", tsv)

# A zero proportion has no place on a log axis. Pin it just below the smallest real
# observation and mark it, rather than dropping the organism silently.
floor_val <- min(d$real_prop[d$real_prop > 0], d$sim_prop[d$sim_prop > 0]) / 3
# TODO (requested 2026-11-25, deferred — do this next time these figures are redrawn):
# abbreviate the genus in point labels, so "Escherichia coli" reads "E. coli". Standard
# binomial convention, and it buys a lot of horizontal room on a crowded log-log panel
# where ggrepel is already fighting for space.
#
#   label = sub("^([A-Za-z])[a-z]+ ", "\\1. ", gsub("_", " ", organism))
#
# Mind the pooled E. coli group: it arrives as "Escherichia_coli" with no strain suffix,
# so the same substitution handles it. Names with no space (a bare genus) are left alone
# by the pattern, which is the behaviour we want.
d <- d %>% mutate(
  censored  = real_prop <= 0 | sim_prop <= 0,
  real_plot = pmax(real_prop, floor_val),
  sim_plot  = pmax(sim_prop,  floor_val),
  fold      = ifelse(censored, NA_real_, sim_prop / real_prop),
  label     = gsub("_", " ", organism)
)

rho <- d$spearman_rho[1]; ccc <- d$ccc_log10[1]; pr <- d$pearson_r_log10[1]
sub <- sprintf("Spearman rho = %s   Pearson r (log10) = %s   Lin's CCC (log10) = %s   n = %d",
               rho, pr, ccc, nrow(d))

lim <- range(c(d$real_plot, d$sim_plot))
lim <- c(lim[1] / 2.5, lim[2] * 2.5)

# On log-log axes a constant fold-change is a line parallel to the diagonal, offset by
# log10(fold). Drawn as guides rather than a shaded polygon, which is illegible at the
# alpha needed to keep it from competing with the points.
p <- ggplot(d, aes(real_plot, sim_plot)) +
  geom_abline(slope = 1, intercept = c(-log10(2), log10(2)),
              linetype = "dotted", colour = "#2980B9", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(fill = log2(fold), shape = censored), size = 3.4, stroke = 0.6,
             colour = "grey25") +
  ggrepel::geom_text_repel(aes(label = label), size = 2.9, colour = "grey20",
                           segment.colour = "grey70", max.overlaps = 20,
                           min.segment.length = 0.1, seed = 42) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24),
                     labels = c(`FALSE` = "observed", `TRUE` = "zero on one side"),
                     name = NULL) +
  scale_fill_gradient2(low = "#C0392B", mid = "#F7F7F7", high = "#1A5276",
                       midpoint = 0, na.value = "grey85",
                       name = expression(log[2]~frac(sim, real))) +
  scale_x_log10(labels = label_percent(accuracy = 0.001), limits = lim) +
  scale_y_log10(labels = label_percent(accuracy = 0.001), limits = lim) +
  annotation_logticks(sides = "bl", colour = "grey60") +
  coord_fixed() +
  labs(title = title, subtitle = sub,
       x = "Real proportion of mapped reads (PacBio CCS)",
       y = "Simulated proportion of mapped reads",
       caption = paste0("Dashed = perfect agreement; dotted = 2-fold either side. The five ",
                        "near-identical E. coli strains are\ncollapsed, since read assignment ",
                        "among them is close to arbitrary for both real and simulated data.")) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 8.5),
        plot.caption  = element_text(colour = "grey45", size = 7.5, hjust = 0),
        legend.key.height = unit(0.8, "lines"))

ggsave(out, p, width = 7.6, height = 7.2, dpi = 200)
cat("wrote", out, "\n")
