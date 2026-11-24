#!/usr/bin/env Rscript
# Scatter of predicted vs measured per-organism ENRICHMENT FACTOR for the D6331 mock
# community — the non-circular validation target (improvements/parameterisation.md B0).
#
#   usage: plot_enrichment.R <enrichment.tsv> <output.png> [title]
#
# The TSV comes from `enrichment.py --tsv`. Enrichment is (share in the enriched arm) /
# (share in the non-enriched arm), so 1.0 is the no-effect line and the quadrants mean
# something: above 1 the panel enriched the organism, below 1 it depleted it. Both axes
# are log10 because the measured values span 63x.
#
# This differs from plot.R, which scatters absolute composition. Composition is
# dominated by how well the input was guessed; enrichment divides that out, so this plot
# shows the capture model and nothing else.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot_enrichment.R <enrichment.tsv> <output.png> [title]")
tsv <- args[1]; out <- args[2]
title <- if (length(args) >= 3) args[3] else
  "D6331 mock community: predicted vs measured enrichment factor"

d <- read.delim(tsv, stringsAsFactors = FALSE)
if (!nrow(d)) stop("no rows in ", tsv)
d$scored <- as.logical(d$scored)

# Organisms excluded for cause are still drawn, hollow, so the reader can see what was
# left out and why rather than having to take the n on trust.
d <- d %>% mutate(
  fold  = ifelse(predicted > 0 & measured > 0, predicted / measured, NA_real_),
  label = gsub("_", " ", organism)
)

sc <- d[d$scored, ]
if (!nrow(sc)) stop("no scored organisms in ", tsv)

# Metrics are recomputed here rather than carried in the TSV so the figure cannot drift
# out of step with the points it draws.
lp <- log10(sc$predicted); lm_ <- log10(sc$measured)
n <- nrow(sc)
rho  <- suppressWarnings(cor(rank(lp), rank(lm_)))
pear <- suppressWarnings(cor(lp, lm_))
rmse <- sqrt(mean((lp - lm_)^2))
cen  <- sqrt(mean(((lp - mean(lp)) - (lm_ - mean(lm_)))^2))
sub <- sprintf(paste0("Spearman rho = %.3f   Pearson r (log10) = %.3f   ",
                      "log10 RMSE = %.3f (%.2fx)   centred = %.3f   n = %d"),
               rho, pear, rmse, 10^rmse, cen, n)

lim <- range(c(d$predicted[d$predicted > 0], d$measured[d$measured > 0]))
lim <- c(lim[1] / 2.5, lim[2] * 2.5)

p <- ggplot(d, aes(measured, predicted)) +
  # 1.0 on either axis separates enrichment from depletion, which is the physically
  # meaningful threshold here — more so than any fold-change band.
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "#B9770E") +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "#B9770E") +
  geom_abline(slope = 1, intercept = c(-log10(2), log10(2)),
              linetype = "dotted", colour = "#2980B9", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(fill = log2(fold), shape = scored, size = input_share),
             stroke = 0.6, colour = "grey25") +
  ggrepel::geom_text_repel(aes(label = label), size = 2.9, colour = "grey20",
                           segment.colour = "grey70", max.overlaps = 20,
                           min.segment.length = 0.1, seed = 42) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24),
                     labels = c(`TRUE` = "scored", `FALSE` = "excluded for cause"),
                     name = NULL) +
  scale_size_continuous(range = c(2, 7), trans = "log10",
                        labels = label_percent(accuracy = 0.01),
                        name = "input share") +
  scale_fill_gradient2(low = "#C0392B", mid = "#F7F7F7", high = "#1A5276",
                       midpoint = 0, na.value = "grey85",
                       name = expression(log[2]~frac(predicted, measured))) +
  scale_x_log10(limits = lim) + scale_y_log10(limits = lim) +
  annotation_logticks(sides = "bl", colour = "grey60") +
  coord_fixed() +
  guides(shape = guide_legend(order = 1), size = guide_legend(order = 2),
         fill = guide_colourbar(order = 3)) +
  labs(title = title, subtitle = sub,
       x = "Measured enrichment (enriched / non-enriched arm)",
       y = "Predicted enrichment (simulated / non-enriched arm)",
       caption = paste0("Input composition comes from the pooled non-enriched arm, so ",
                        "input and target are different libraries.\nOrange lines mark ",
                        "enrichment = 1 (no effect); dashed = perfect agreement; ",
                        "dotted = 2-fold either side.")) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey30", size = 8.5),
        plot.caption  = element_text(colour = "grey45", size = 7.5, hjust = 0),
        legend.key.height = unit(0.8, "lines"))

ggsave(out, p, width = 8.0, height = 7.4, dpi = 200)
cat("wrote", out, "\n")
