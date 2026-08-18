#!/usr/bin/env Rscript
# BaitBench Diagrams
# Generates diagrams for prepare modes and pipeline flowcharts.
#
# Usage: Rscript docs/prepare_diagrams.R [output_dir] [diagrams]
#
# Arguments:
#   output_dir  Directory for PNGs (default: ".")
#   diagrams    Which diagrams to render (default: all). Accepts:
#               - Single number: 5
#               - Comma-separated: 1,3,5
#               - Ranges: 2-4
#               - Mixed: 1,3,5-6
#
# Diagrams:
#   1  prepare_mode1_standard_nosample.png
#   2  prepare_mode2_standard_sample.png
#   3  prepare_mode3_genomes_nosample.png
#   4  prepare_mode4_genomes_sample.png
#   5  pipeline_overview.png
#   6  pipeline_detailed.png
#   7  paper_workflow_overview.png        — full tool workflow (design→assess→simulate→report)
#   8  paper_thermodynamic_scoring.png   — TNN ΔG scoring illustration
#   9  paper_fragment_sampling.png       — two-level multinomial sampling schematic
#  10  paper_build_probes.png            — build-probes pipeline steps
#  11  paper_assess_probes.png           — assess-probes pipeline steps
#  12  paper_full_simulation_pipeline.png — genome-mode simulation pipeline, full detail
#                                       (prepare → probe alignment → thermodynamic
#                                       scoring → fragment sampling → sequence
#                                       simulator options → filter → map/list →
#                                       metrics → report). Intentionally dense —
#                                       not meant to be readable at a glance.

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Parse diagram selection (e.g. "5", "2-4", "1,3,5-6")
parse_selection <- function(spec, total = 12) {
  if (is.na(spec) || spec == "" || spec == "all") return(seq_len(total))
  parts <- unlist(strsplit(spec, ","))
  nums <- integer(0)
  for (p in parts) {
    p <- trimws(p)
    if (grepl("^\\d+-\\d+$", p)) {
      bounds <- as.integer(unlist(strsplit(p, "-")))
      nums <- c(nums, seq(bounds[1], bounds[2]))
    } else if (grepl("^\\d+$", p)) {
      nums <- c(nums, as.integer(p))
    } else {
      stop("Invalid diagram specifier: ", p)
    }
  }
  sort(unique(nums))
}

diagram_spec <- if (length(args) >= 2) args[2] else "all"
selected <- parse_selection(diagram_spec)
message("Rendering diagrams: ", paste(selected, collapse = ", "))

# Poster mode: `poster` as the 3rd argument emits high-resolution,
# transparent-background copies prefixed with poster_ instead of the
# normal web-sized renders.
poster_mode  <- length(args) >= 3 && tolower(args[3]) %in% c("poster", "true", "yes")
POSTER_SCALE <- 3          # multiplier applied to each diagram base width
POSTER_MAX_PX <- 140e6     # pixel-count ceiling, keeps rsvg inside memory

if (poster_mode) message("Poster mode: high-res, transparent background, poster_ prefix")

# Re-render a graph with a transparent canvas.
# The DOT source lives in graph$x$diagram, so swap the graph-level bgcolor
# there rather than post-processing the SVG. Node `fillcolor='white'` and HTML
# table `BGCOLOR=` attributes do not contain the substring "bgcolor='white'",
# so intentional white fills survive untouched.
make_transparent <- function(graph) {
  dot <- graph$x$diagram
  dot <- gsub("bgcolor='white'", "bgcolor='transparent'", dot, fixed = TRUE)
  dot <- gsub('bgcolor="white"', 'bgcolor="transparent"', dot, fixed = TRUE)
  grViz(dot)
}

# Read the intrinsic point dimensions graphviz wrote into the SVG header so the
# pixel ceiling can be applied without guessing the aspect ratio.
svg_aspect <- function(svg) {
  w <- as.numeric(sub('.*<svg[^>]*width="([0-9.]+)pt".*', '\\1', svg))
  h <- as.numeric(sub('.*<svg[^>]*height="([0-9.]+)pt".*', '\\1', svg))
  if (is.na(w) || is.na(h) || w <= 0) return(NA_real_)
  h / w
}

save_poster <- function(graph, filename, width) {
  svg <- export_svg(make_transparent(graph))
  target <- width * POSTER_SCALE
  aspect <- svg_aspect(svg)
  if (!is.na(aspect)) {
    px <- target * (target * aspect)
    if (px > POSTER_MAX_PX) {
      target <- floor(sqrt(POSTER_MAX_PX / aspect))
      message("  (capped to ", target, " px wide to stay within memory)")
    }
  }
  out <- file.path(outdir, paste0("poster_", sub("^paper_", "", filename)))
  rsvg_png(charToRaw(svg), file = out, width = target)
  message("Saved: ", out, "  [", target, " px wide]")
}

save_diagram <- function(graph, filename, width = 3200) {
  if (poster_mode) {
    save_poster(graph, filename, width)
    return(invisible(NULL))
  }
  svg <- export_svg(graph)
  rsvg_png(charToRaw(svg), file = file.path(outdir, filename), width = width)
  message("Saved: ", file.path(outdir, filename))
}


# ============================================================================
# Diagram 1: Standard Mode, No --sample
# ============================================================================
if (1 %in% selected) {
diagram1 <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>Mode 1: Standard (no --genomes, no --sample)</B><BR/>All targets are in sample with weight 1.0>,
         labelloc=t, fontsize=20, pad=0.5, nodesep=0.5, ranksep=1.0]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box, margin='0.15,0.08']
  edge [fontname='Helvetica', fontsize=10]

  // --- INPUT FILES ---
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=14; fontcolor='#7F8C8D'

    targets [label=<<B>targets.fa</B><BR/><FONT POINT-SIZE='10'>Target sequences (e.g. viruses)<BR/>IDs: T1, T2, T3</FONT>>,
             fillcolor='#D6EAF8', color='#4A90D9', penwidth=2]

    distractors [label=<<B>distractors.fa</B><BR/><FONT POINT-SIZE='10'>Background sequences<BR/>IDs: D1, D2, ... Dn</FONT>>,
                 fillcolor='#FADBD8', color='#E74C3C', penwidth=2]

    params [label=<<B>Parameters</B><BR/><FONT POINT-SIZE='10'>--distractor-fraction 0.9<BR/>(or --ct 25)</FONT>>,
            fillcolor='#F2F3F4', color='#95A5A6', shape=note, penwidth=1]
  }

  // --- PREPARE PROCESS ---
  prepare [label=<<B>baitbench prepare</B><BR/><FONT POINT-SIZE='10'>Combine references<BR/>Generate weights</FONT>>,
           fillcolor='#34495E', fontcolor=white, shape=box, penwidth=2]

  // --- OUTPUTS ---
  subgraph cluster_outputs {
    label=<<B>Prepare Outputs</B>>; style=dashed; color='#27AE60'; fontsize=14; fontcolor='#27AE60'

    combined [label=<<B>combined_reference.fa</B><BR/><FONT POINT-SIZE='10'>targets.fa + distractors.fa<BR/>Contains: T1, T2, T3, D1, D2, ... Dn</FONT>>,
              fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    weights [label=<<B>weights.txt</B><BR/><FONT POINT-SIZE='10'>T1  1.0   (sample, default)<BR/>T2  1.0   (sample, default)<BR/>T3  1.0   (sample, default)<BR/>D1  0.009 (distractor weight)<BR/>D2  0.009 (distractor weight)<BR/>...</FONT>>,
             fillcolor='#D5F5E3', color='#2ECC71', penwidth=2]

    idlists [label=<<B>ID Lists</B><BR/><FONT POINT-SIZE='10'>targets.txt: T1, T2, T3<BR/>distractors.txt: D1, D2, ...Dn<BR/>sample.txt: T1, T2, T3  (= all targets)</FONT>>,
             fillcolor='#D1F2EB', color='#1ABC9C', penwidth=2]
  }

  // --- SIMULATE (right side) ---
  simulate [label=<<B>baitbench simulate</B><BR/><FONT POINT-SIZE='10'>Weighted random fragments</FONT>>,
            fillcolor='#2C3E50', fontcolor=white, shape=box, penwidth=2]

  // --- FRAGMENT COMPOSITION (right side) ---
  subgraph cluster_fragments {
    label=<<B>Fragment Composition (fragments.fa)</B>>; style=dashed; color='#2980B9'; fontsize=14; fontcolor='#2980B9'

    frag_comp [label=<
      <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0' CELLPADDING='6'>
        <TR><TD BGCOLOR='#D6EAF8' COLSPAN='2'><B>~10% from targets (sample)</B></TD></TR>
        <TR><TD BGCOLOR='#D6EAF8'>T1 fragments</TD><TD BGCOLOR='#D6EAF8'>weight 1.0</TD></TR>
        <TR><TD BGCOLOR='#D6EAF8'>T2 fragments</TD><TD BGCOLOR='#D6EAF8'>weight 1.0</TD></TR>
        <TR><TD BGCOLOR='#D6EAF8'>T3 fragments</TD><TD BGCOLOR='#D6EAF8'>weight 1.0</TD></TR>
        <TR><TD BGCOLOR='#FADBD8' COLSPAN='2'><B>~90% from distractors</B></TD></TR>
        <TR><TD BGCOLOR='#FADBD8'>D1..Dn fragments</TD><TD BGCOLOR='#FADBD8'>weight calculated<BR/>per distractor fraction</TD></TR>
      </TABLE>>,
      shape=plaintext]
  }

  // --- EDGES ---
  targets -> prepare
  distractors -> prepare
  params -> prepare [style=dashed]

  prepare -> combined
  prepare -> weights
  prepare -> idlists

  combined -> simulate
  weights -> simulate

  simulate -> frag_comp
}
")

save_diagram(diagram1, "prepare_mode1_standard_nosample.png")
}


# ============================================================================
# Diagram 2: Standard Mode, With --sample
# ============================================================================
if (2 %in% selected) {
diagram2 <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>Mode 2: Standard (no --genomes, with --sample)</B><BR/>Only sample targets generate fragments; non-sample targets get weight 0>,
         labelloc=t, fontsize=20, pad=0.5, nodesep=0.5, ranksep=1.0]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box, margin='0.15,0.08']
  edge [fontname='Helvetica', fontsize=10]

  // --- INPUT FILES ---
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=14; fontcolor='#7F8C8D'

    targets [label=<<B>targets.fa</B><BR/><FONT POINT-SIZE='10'>Target sequences<BR/>IDs: T1, T2, T3</FONT>>,
             fillcolor='#D6EAF8', color='#4A90D9', penwidth=2]

    distractors [label=<<B>distractors.fa</B><BR/><FONT POINT-SIZE='10'>Background sequences<BR/>IDs: D1, D2, ... Dn</FONT>>,
                 fillcolor='#FADBD8', color='#E74C3C', penwidth=2]

    sample [label=<<B>--sample manifest.tsv</B><BR/><FONT POINT-SIZE='10'>T1  1.0<BR/>T2  5.0<BR/>(T3 is NOT in sample)</FONT>>,
            fillcolor='#FDEBD0', color='#F39C12', penwidth=2, shape=note]

    params [label=<<B>Parameters</B><BR/><FONT POINT-SIZE='10'>--distractor-fraction 0.9</FONT>>,
            fillcolor='#F2F3F4', color='#95A5A6', shape=note, penwidth=1]
  }

  // --- PREPARE PROCESS ---
  prepare [label=<<B>baitbench prepare</B><BR/><FONT POINT-SIZE='10'>Validate sample IDs in targets<BR/>Combine references<BR/>Generate weights</FONT>>,
           fillcolor='#34495E', fontcolor=white, shape=box, penwidth=2]

  // --- OUTPUTS ---
  subgraph cluster_outputs {
    label=<<B>Prepare Outputs</B>>; style=dashed; color='#27AE60'; fontsize=14; fontcolor='#27AE60'

    combined [label=<<B>combined_reference.fa</B><BR/><FONT POINT-SIZE='10'>targets.fa + distractors.fa<BR/>Contains: T1, T2, T3, D1, D2, ... Dn</FONT>>,
              fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    weights [label=<<B>weights.txt</B><BR/><FONT POINT-SIZE='10'>T1  1.0    (sample)<BR/>T2  5.0    (sample, high weight)<BR/><FONT COLOR='#95A5A6'>T3  0.0    (non-sample, no fragments!)</FONT><BR/>D1  0.054  (distractor weight)<BR/>D2  0.054  (distractor weight)<BR/>...</FONT>>,
             fillcolor='#D5F5E3', color='#2ECC71', penwidth=2]

    idlists [label=<<B>ID Lists</B><BR/><FONT POINT-SIZE='10'>targets.txt: T1, T2, T3<BR/>distractors.txt: D1, D2, ...Dn<BR/>sample.txt: T1, T2  (subset!)</FONT>>,
             fillcolor='#D1F2EB', color='#1ABC9C', penwidth=2]
  }

  // --- SIMULATE (right side) ---
  simulate [label=<<B>baitbench simulate</B><BR/><FONT POINT-SIZE='10'>Weighted random fragments</FONT>>,
            fillcolor='#2C3E50', fontcolor=white, shape=box, penwidth=2]

  // --- FRAGMENT COMPOSITION (right side) ---
  subgraph cluster_fragments {
    label=<<B>Fragment Composition (fragments.fa)</B>>; style=dashed; color='#2980B9'; fontsize=14; fontcolor='#2980B9'

    frag_comp [label=<
      <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0' CELLPADDING='6'>
        <TR><TD BGCOLOR='#D6EAF8' COLSPAN='2'><B>~10% from sample targets</B></TD></TR>
        <TR><TD BGCOLOR='#D6EAF8'>T1 fragments</TD><TD BGCOLOR='#D6EAF8'>weight 1.0 (~1/6 of target share)</TD></TR>
        <TR><TD BGCOLOR='#D6EAF8'>T2 fragments</TD><TD BGCOLOR='#D6EAF8'>weight 5.0 (~5/6 of target share)</TD></TR>
        <TR><TD BGCOLOR='#E8E8E8' COLSPAN='2'><B>0% from non-sample targets</B></TD></TR>
        <TR><TD BGCOLOR='#E8E8E8'><FONT COLOR='#95A5A6'>T3 fragments</FONT></TD><TD BGCOLOR='#E8E8E8'><FONT COLOR='#95A5A6'>weight 0.0 (none generated)</FONT></TD></TR>
        <TR><TD BGCOLOR='#FADBD8' COLSPAN='2'><B>~90% from distractors</B></TD></TR>
        <TR><TD BGCOLOR='#FADBD8'>D1..Dn fragments</TD><TD BGCOLOR='#FADBD8'>weight per distractor fraction</TD></TR>
      </TABLE>>,
      shape=plaintext]
  }

  // --- EDGES ---
  targets -> prepare
  distractors -> prepare
  sample -> prepare
  params -> prepare [style=dashed]

  prepare -> combined
  prepare -> weights
  prepare -> idlists

  combined -> simulate
  weights -> simulate

  simulate -> frag_comp
}
")

save_diagram(diagram2, "prepare_mode2_standard_sample.png")
}


# ============================================================================
# Diagram 3: Genome Mode, No --sample
# ============================================================================
if (3 %in% selected) {
diagram3 <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>Mode 3: Genome Mode (with --genomes, no --sample)</B><BR/>Fragments come from full genomes; reads map back to targets>,
         labelloc=t, fontsize=20, pad=0.5, nodesep=0.5, ranksep=1.0]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box, margin='0.15,0.08']
  edge [fontname='Helvetica', fontsize=10]

  // --- INPUT FILES ---
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=14; fontcolor='#7F8C8D'

    genomes [label=<<B>genomes.fa</B><BR/><FONT POINT-SIZE='10'>Full genome sequences<BR/>(e.g. whole bacteria)<BR/>IDs: G1, G2</FONT>>,
             fillcolor='#E8DAEF', color='#8E44AD', penwidth=2]

    targets [label=<<B>targets.fa</B><BR/><FONT POINT-SIZE='10'>Probe target subsequences<BR/>(e.g. 16S gene from each)<BR/>IDs: G1|16S, G2|ompB</FONT>>,
             fillcolor='#D6EAF8', color='#4A90D9', penwidth=2]

    distractors [label=<<B>distractors.fa</B><BR/><FONT POINT-SIZE='10'>Background sequences<BR/>IDs: D1, D2, ... Dn</FONT>>,
                 fillcolor='#FADBD8', color='#E74C3C', penwidth=2]

    params [label=<<B>Parameters</B><BR/><FONT POINT-SIZE='10'>--distractor-fraction 0.9</FONT>>,
            fillcolor='#F2F3F4', color='#95A5A6', shape=note, penwidth=1]
  }

  // --- PREPARE PROCESS ---
  prepare [label=<<B>baitbench prepare</B><BR/><FONT POINT-SIZE='10'>Build two references<BR/>Auto-link genomes to targets<BR/>All genomes = sample (wt 1.0)</FONT>>,
           fillcolor='#34495E', fontcolor=white, shape=box, penwidth=2]

  // --- OUTPUTS ---
  subgraph cluster_outputs {
    label=<<B>Prepare Outputs</B>>; style=dashed; color='#27AE60'; fontsize=14; fontcolor='#27AE60'

    combined [label=<<B>combined_reference.fa</B><BR/><FONT POINT-SIZE='10'><B>genomes.fa</B> + distractors.fa<BR/>For fragment generation<BR/>Contains: G1, G2, D1, ... Dn</FONT>>,
              fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    mapping [label=<<B>mapping_reference.fa</B><BR/><FONT POINT-SIZE='10'><B>targets.fa</B> + distractors.fa<BR/>For read mapping (later)<BR/>Contains: G1|16S, G2|ompB, D1, ... Dn</FONT>>,
             fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    weights [label=<<B>weights.txt</B><BR/><FONT POINT-SIZE='10'>G1  1.0   (sample, default)<BR/>G2  1.0   (sample, default)<BR/>D1  0.009 (distractor weight)<BR/>D2  0.009 (distractor weight)<BR/>...</FONT>>,
             fillcolor='#D5F5E3', color='#2ECC71', penwidth=2]

    stmap [label=<<B>sample_target_map.txt</B><BR/><FONT POINT-SIZE='10'>Auto-linked by prefix match:<BR/>G1 &rarr; G1|16S<BR/>G2 &rarr; G2|ompB</FONT>>,
           fillcolor='#FDEBD0', color='#E67E22', penwidth=2]

    idlists [label=<<B>ID Lists</B><BR/><FONT POINT-SIZE='10'>genomes.txt: G1, G2<BR/>targets.txt: G1|16S, G2|ompB<BR/>distractors.txt: D1, D2, ...Dn<BR/>sample.txt: G1, G2  (= all genomes)</FONT>>,
             fillcolor='#D1F2EB', color='#1ABC9C', penwidth=2]
  }

  // --- SIMULATE (right side) ---
  simulate [label=<<B>baitbench simulate</B><BR/><FONT POINT-SIZE='10'>Weighted random fragments<BR/>from <B>full genomes</B></FONT>>,
            fillcolor='#2C3E50', fontcolor=white, shape=box, penwidth=2]

  // --- FRAGMENT COMPOSITION (right side) ---
  subgraph cluster_fragments {
    label=<<B>Fragment Composition (fragments.fa)</B>>; style=dashed; color='#2980B9'; fontsize=14; fontcolor='#2980B9'

    frag_comp [label=<
      <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0' CELLPADDING='6'>
        <TR><TD BGCOLOR='#E8DAEF' COLSPAN='2'><B>~10% from genomes (sample)</B></TD></TR>
        <TR><TD BGCOLOR='#E8DAEF'>G1 fragments</TD><TD BGCOLOR='#E8DAEF'>weight 1.0 (full genome!<BR/>most frags NOT on target)</TD></TR>
        <TR><TD BGCOLOR='#E8DAEF'>G2 fragments</TD><TD BGCOLOR='#E8DAEF'>weight 1.0 (full genome!<BR/>most frags NOT on target)</TD></TR>
        <TR><TD BGCOLOR='#FADBD8' COLSPAN='2'><B>~90% from distractors</B></TD></TR>
        <TR><TD BGCOLOR='#FADBD8'>D1..Dn fragments</TD><TD BGCOLOR='#FADBD8'>weight per distractor fraction</TD></TR>
      </TABLE>>,
      shape=plaintext]

    note [label=<<FONT POINT-SIZE='10'><I>Fragments span the full genome, not just<BR/>probe target regions. Only a small fraction<BR/>of genome fragments overlap the target<BR/>gene (e.g. 16S). Probes must capture<BR/>these from off-target genomic fragments.</I></FONT>>,
          shape=note, fillcolor='#FEF9E7', color='#F1C40F', penwidth=1]
  }

  // --- EDGES ---
  genomes -> prepare
  targets -> prepare
  distractors -> prepare
  params -> prepare [style=dashed]

  prepare -> combined
  prepare -> mapping
  prepare -> weights
  prepare -> stmap
  prepare -> idlists

  combined -> simulate
  weights -> simulate

  simulate -> frag_comp
}
")

save_diagram(diagram3, "prepare_mode3_genomes_nosample.png")
}


# ============================================================================
# Diagram 4: Genome Mode, With --sample
# ============================================================================
if (4 %in% selected) {
diagram4 <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>Mode 4: Genome Mode (with --genomes and --sample)</B><BR/>Only sample genomes generate fragments; explicit or auto-linked target mapping>,
         labelloc=t, fontsize=20, pad=0.5, nodesep=0.5, ranksep=1.0]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box, margin='0.15,0.08']
  edge [fontname='Helvetica', fontsize=10]

  // --- INPUT FILES ---
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=14; fontcolor='#7F8C8D'

    genomes [label=<<B>genomes.fa</B><BR/><FONT POINT-SIZE='10'>Full genome sequences<BR/>IDs: G1, G2, G3</FONT>>,
             fillcolor='#E8DAEF', color='#8E44AD', penwidth=2]

    targets [label=<<B>targets.fa</B><BR/><FONT POINT-SIZE='10'>Probe target subsequences<BR/>IDs: G1|16S, G2|ompB, G3|gltA</FONT>>,
             fillcolor='#D6EAF8', color='#4A90D9', penwidth=2]

    distractors [label=<<B>distractors.fa</B><BR/><FONT POINT-SIZE='10'>Background sequences<BR/>IDs: D1, D2, ... Dn</FONT>>,
                 fillcolor='#FADBD8', color='#E74C3C', penwidth=2]

    sample [label=<<B>--sample manifest.tsv</B><BR/><FONT POINT-SIZE='10'>G1  1.0<BR/>G2  3.0<BR/>(G3 is NOT in sample)</FONT>>,
            fillcolor='#FDEBD0', color='#F39C12', penwidth=2, shape=note]

    stmap_in [label=<<B>--sample-target-map</B><BR/><FONT POINT-SIZE='10'>(optional explicit TSV)<BR/>G1  G1|16S<BR/>G2  G2|ompB<BR/>G3  G3|gltA<BR/><I>Or: auto-linked by prefix</I></FONT>>,
              fillcolor='#FDEBD0', color='#E67E22', penwidth=2, shape=note]

    params [label=<<B>Parameters</B><BR/><FONT POINT-SIZE='10'>--distractor-fraction 0.9</FONT>>,
            fillcolor='#F2F3F4', color='#95A5A6', shape=note, penwidth=1]
  }

  // --- PREPARE PROCESS ---
  prepare [label=<<B>baitbench prepare</B><BR/><FONT POINT-SIZE='10'>Validate sample IDs in genomes<BR/>Resolve genome &rarr; target maps<BR/>Build two references<BR/>Generate weights</FONT>>,
           fillcolor='#34495E', fontcolor=white, shape=box, penwidth=2]

  // --- OUTPUTS ---
  subgraph cluster_outputs {
    label=<<B>Prepare Outputs</B>>; style=dashed; color='#27AE60'; fontsize=14; fontcolor='#27AE60'

    combined [label=<<B>combined_reference.fa</B><BR/><FONT POINT-SIZE='10'><B>genomes.fa</B> + distractors.fa<BR/>For fragment generation<BR/>Contains: G1, G2, G3, D1, ... Dn</FONT>>,
              fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    mapping [label=<<B>mapping_reference.fa</B><BR/><FONT POINT-SIZE='10'><B>targets.fa</B> + distractors.fa<BR/>For read mapping (later)<BR/>Contains: G1|16S, G2|ompB, G3|gltA,<BR/>D1, ... Dn</FONT>>,
             fillcolor='#D5F5E3', color='#27AE60', penwidth=2]

    weights [label=<<B>weights.txt</B><BR/><FONT POINT-SIZE='10'>G1  1.0   (sample)<BR/>G2  3.0   (sample, high weight)<BR/><FONT COLOR='#95A5A6'>G3  0.0   (non-sample, no frags!)</FONT><BR/>D1  0.036 (distractor weight)<BR/>D2  0.036 (distractor weight)<BR/>...</FONT>>,
             fillcolor='#D5F5E3', color='#2ECC71', penwidth=2]

    stmap_out [label=<<B>sample_target_map.txt</B><BR/><FONT POINT-SIZE='10'>G1 &rarr; G1|16S<BR/>G2 &rarr; G2|ompB<BR/>G3 &rarr; G3|gltA<BR/>(used by metrics step)</FONT>>,
               fillcolor='#FDEBD0', color='#E67E22', penwidth=2]

    idlists [label=<<B>ID Lists</B><BR/><FONT POINT-SIZE='10'>genomes.txt: G1, G2, G3<BR/>targets.txt: G1|16S, G2|ompB, G3|gltA<BR/>distractors.txt: D1, D2, ...Dn<BR/>sample.txt: G1, G2  (subset!)</FONT>>,
             fillcolor='#D1F2EB', color='#1ABC9C', penwidth=2]
  }

  // --- SIMULATE (right side) ---
  simulate [label=<<B>baitbench simulate</B><BR/><FONT POINT-SIZE='10'>Weighted random fragments<BR/>from <B>full genomes</B></FONT>>,
            fillcolor='#2C3E50', fontcolor=white, shape=box, penwidth=2]

  // --- FRAGMENT COMPOSITION (right side) ---
  subgraph cluster_fragments {
    label=<<B>Fragment Composition (fragments.fa)</B>>; style=dashed; color='#2980B9'; fontsize=14; fontcolor='#2980B9'

    frag_comp [label=<
      <TABLE BORDER='0' CELLBORDER='1' CELLSPACING='0' CELLPADDING='6'>
        <TR><TD BGCOLOR='#E8DAEF' COLSPAN='2'><B>~10% from sample genomes</B></TD></TR>
        <TR><TD BGCOLOR='#E8DAEF'>G1 fragments</TD><TD BGCOLOR='#E8DAEF'>weight 1.0 (~1/4 of target share)<BR/>full genome fragments</TD></TR>
        <TR><TD BGCOLOR='#E8DAEF'>G2 fragments</TD><TD BGCOLOR='#E8DAEF'>weight 3.0 (~3/4 of target share)<BR/>full genome fragments</TD></TR>
        <TR><TD BGCOLOR='#E8E8E8' COLSPAN='2'><B>0% from non-sample genomes</B></TD></TR>
        <TR><TD BGCOLOR='#E8E8E8'><FONT COLOR='#95A5A6'>G3 fragments</FONT></TD><TD BGCOLOR='#E8E8E8'><FONT COLOR='#95A5A6'>weight 0.0 (none generated)</FONT></TD></TR>
        <TR><TD BGCOLOR='#FADBD8' COLSPAN='2'><B>~90% from distractors</B></TD></TR>
        <TR><TD BGCOLOR='#FADBD8'>D1..Dn fragments</TD><TD BGCOLOR='#FADBD8'>weight per distractor fraction</TD></TR>
      </TABLE>>,
      shape=plaintext]

    note [label=<<FONT POINT-SIZE='10'><I>G3 is in genomes.fa but NOT in --sample,<BR/>so weight=0, no fragments generated.<BR/>G3|gltA still exists in mapping_reference.fa;<BR/>reads from other sources could mis-map<BR/>to it (counted as FP_target by metrics).</I></FONT>>,
          shape=note, fillcolor='#FEF9E7', color='#F1C40F', penwidth=1]
  }

  // --- EDGES ---
  genomes -> prepare
  targets -> prepare
  distractors -> prepare
  sample -> prepare
  stmap_in -> prepare [style=dashed, label=<<FONT POINT-SIZE='9'>optional</FONT>>]
  params -> prepare [style=dashed]

  prepare -> combined
  prepare -> mapping
  prepare -> weights
  prepare -> stmap_out
  prepare -> idlists

  combined -> simulate
  weights -> simulate

  simulate -> frag_comp
}
")

save_diagram(diagram4, "prepare_mode4_genomes_sample.png")
}


# ============================================================================
# Diagram 5: Pipeline High-Level Overview
# Clean horizontal flow for landscape presentation slides.
# ============================================================================
if (5 %in% selected) {
pipeline_overview <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>BaitBench Pipeline Overview</B>>,
         labelloc=t, fontsize=24, pad='0.5,0.3', nodesep=0.4, ranksep=0.7]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box,
        margin='0.15,0.08', penwidth=1.5]
  edge [fontname='Helvetica', fontsize=9, color='#5D6D7E', penwidth=1.5]

  // --- USER INPUTS (top-left) ---
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=rounded; color='#BDC3C7'; fontsize=13; fontcolor='#7F8C8D'
    margin=15

    targets [label=<<B>Targets</B><BR/><FONT POINT-SIZE='9'>Target sequences</FONT>>,
             fillcolor='#D6EAF8', color='#4A90D9']
    distractors [label=<<B>Distractors</B><BR/><FONT POINT-SIZE='9'>Background sequences</FONT>>,
                 fillcolor='#FADBD8', color='#E74C3C']
    probes [label=<<B>Probes</B><BR/><FONT POINT-SIZE='9'>Capture probes</FONT>>,
            fillcolor='#FCF3CF', color='#D4AC0D']
    sample_in [label=<<B>Sample</B><BR/><FONT POINT-SIZE='9'>Manifest (optional)</FONT>>,
               fillcolor='#FDEBD0', color='#F39C12', style='filled,dashed']
    host_in [label=<<B>Host</B><BR/><FONT POINT-SIZE='9'>Host genome (optional)</FONT>>,
             fillcolor='#F2F3F4', color='#95A5A6', style='filled,dashed']
  }

  // --- PIPELINE STEPS: main chain ---
  prepare [label=<<B>1. Prepare</B><BR/><FONT POINT-SIZE='9'>Combine references,<BR/>generate weights</FONT>>,
           fillcolor='#34495E', fontcolor='white', shape=box]

  simulate [label=<<B>2. Simulate</B><BR/><FONT POINT-SIZE='9'>Weighted random<BR/>fragments</FONT>>,
            fillcolor='#2C3E50', fontcolor='white']

  capture [label=<<B>3. Capture</B><BR/><FONT POINT-SIZE='9'>Probe hybridization<BR/>(minimap2 / BLAST)</FONT>>,
           fillcolor='#1A5276', fontcolor='white']

  enrich [label=<<B>4. Enrich</B><BR/><FONT POINT-SIZE='9'>Fold enrichment</FONT>>,
          fillcolor='#21618C', fontcolor='white', style='filled,dashed']

  sequence [label=<<B>5. Sequence</B><BR/><FONT POINT-SIZE='9'>Trim to read length,<BR/>subsample</FONT>>,
            fillcolor='#1B4F72', fontcolor='white']

  filter_step [label=<<B>6. Filter</B><BR/><FONT POINT-SIZE='9'>Remove host reads</FONT>>,
               fillcolor='#1B4F72', fontcolor='white', style='filled,dashed']

  map_step [label=<<B>7. Map</B><BR/><FONT POINT-SIZE='9'>Align reads to<BR/>references</FONT>>,
            fillcolor='#154360', fontcolor='white']

  list_step [label=<<B>8. List</B><BR/><FONT POINT-SIZE='9'>Count reads<BR/>per reference</FONT>>,
             fillcolor='#154360', fontcolor='white']

  metrics [label=<<B>9. Metrics</B><BR/><FONT POINT-SIZE='9'>TP / FP / FN / TN<BR/>classification</FONT>>,
           fillcolor='#0B5345', fontcolor='white']

  report [label=<<B>10. Report</B><BR/><FONT POINT-SIZE='9'>HTML report<BR/>(R / ggplot2)</FONT>>,
          fillcolor='#0B5345', fontcolor='white']

  // --- OUTPUT ---
  output [label=<<B>report.html</B><BR/><FONT POINT-SIZE='9'>Interactive HTML with<BR/>figures &amp; tables</FONT>>,
          fillcolor='#D5F5E3', color='#27AE60', penwidth=2.5, shape=note]

  // --- MAIN CHAIN (high-weight edges keep it straight) ---
  prepare -> simulate [weight=10]
  simulate -> capture [weight=10]
  capture -> enrich [weight=10, style=dashed, label=<<FONT POINT-SIZE='8'>optional</FONT>>]
  enrich -> sequence [weight=10]
  sequence -> filter_step [weight=10, style=dashed, label=<<FONT POINT-SIZE='8'>optional</FONT>>]
  filter_step -> map_step [weight=10]
  map_step -> list_step [weight=10]
  list_step -> metrics [weight=10]
  metrics -> report [weight=10]
  report -> output [weight=10]

  // --- INPUT EDGES ---
  targets -> prepare [weight=2]
  distractors -> prepare [weight=2]
  sample_in -> prepare [style=dashed, weight=1]
  probes -> capture [weight=1]
  host_in -> filter_step [style=dashed, weight=1]

  // --- LEGEND (bottom) ---
  subgraph cluster_legend {
    label=''; style=invis; margin=6
    node [fontsize=9, margin='0.06,0.03']

    leg_req [label='Required step', fillcolor='#1A5276', fontcolor='white']
    leg_opt [label='Optional step', fillcolor='#1A5276', fontcolor='white', style='filled,dashed']
    leg_inp [label='User input', fillcolor='#D6EAF8', color='#4A90D9']

    leg_req -> leg_opt -> leg_inp [style=invis]
  }
}
")

save_diagram(pipeline_overview, "pipeline_overview.png")
}


# ============================================================================
# Diagram 6: Pipeline Detailed with All Input/Output Files
# Landscape LR layout: step → file → step → file chain.
# Secondary inputs (probes, host, reference reuse, ID lists) shown as
# lighter edges feeding into steps from above/below.
# ============================================================================
if (6 %in% selected) {
pipeline_detailed <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>BaitBench Pipeline — Detailed Data Flow</B>>,
         labelloc=t, fontsize=24, pad='0.5,0.3', nodesep=0.3, ranksep=0.5]
  node [fontname='Helvetica', fontsize=10, style=filled, margin='0.10,0.05', penwidth=1.5]
  edge [fontname='Helvetica', fontsize=8, color='#5D6D7E', penwidth=1.2]

  // ==================================================================
  // USER INPUTS
  // ==================================================================
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=rounded; color='#BDC3C7'; fontsize=12; fontcolor='#7F8C8D'
    margin=12
    node [shape=note, fontsize=10]

    in_targets [label=<<B>targets.fa</B>>, fillcolor='#D6EAF8', color='#4A90D9']
    in_distractors [label=<<B>distractors.fa</B>>, fillcolor='#FADBD8', color='#E74C3C']
    in_sample [label=<<B>sample.tsv</B><BR/><FONT POINT-SIZE='8'>(optional)</FONT>>,
               fillcolor='#FDEBD0', color='#F39C12']
    in_probes [label=<<B>probes.fa</B>>, fillcolor='#FCF3CF', color='#D4AC0D']
    in_host [label=<<B>host.fa</B><BR/><FONT POINT-SIZE='8'>(optional)</FONT>>,
             fillcolor='#F2F3F4', color='#95A5A6']
  }

  // ==================================================================
  // PIPELINE: step → file → step → file  (main chain)
  // ==================================================================

  // --- PREPARE ---
  prepare [label=<<B>1. Prepare</B><BR/><FONT POINT-SIZE='8'>Combine references<BR/>generate weights</FONT>>,
           shape=box, fillcolor='#34495E', fontcolor='white']

  f_combined [label=<<B>combined_reference.fa</B>>,
              shape=note, fillcolor='#D5F5E3', color='#27AE60']
  f_weights [label=<<B>weights.txt</B>>,
             shape=note, fillcolor='#D5F5E3', color='#27AE60']
  f_idlists [label=<<B>targets.txt<BR/>distractors.txt<BR/>sample.txt</B>>,
             shape=note, fillcolor='#D1F2EB', color='#1ABC9C']

  // --- SIMULATE ---
  simulate [label=<<B>2. Simulate</B><BR/><FONT POINT-SIZE='8'>Weighted random<BR/>fragments</FONT>>,
            shape=box, fillcolor='#2C3E50', fontcolor='white']

  f_fragments [label=<<B>fragments.fa</B>>,
               shape=note, fillcolor='#D5F5E3', color='#27AE60']

  // --- CAPTURE ---
  capture [label=<<B>3. Capture</B><BR/><FONT POINT-SIZE='8'>minimap2 / BLAST</FONT>>,
           shape=box, fillcolor='#1A5276', fontcolor='white']

  f_captured [label=<<B>captured.fa</B>>,
              shape=note, fillcolor='#D5F5E3', color='#27AE60']

  // --- ENRICH (optional) ---
  enrich [label=<<B>4. Enrich</B><BR/><FONT POINT-SIZE='8'>(optional)</FONT>>,
          shape=box, fillcolor='#21618C', fontcolor='white', style='filled,dashed']

  f_enriched [label=<<B>enriched.fa</B>>,
              shape=note, fillcolor='#D5F5E3', color='#2ECC71', style='filled,dashed']

  // --- SEQUENCE ---
  sequence [label=<<B>5. Sequence</B><BR/><FONT POINT-SIZE='8'>Trim &amp; subsample</FONT>>,
            shape=box, fillcolor='#1B4F72', fontcolor='white']

  f_reads [label=<<B>reads.fa</B>>,
           shape=note, fillcolor='#D5F5E3', color='#27AE60']

  // --- FILTER (optional) ---
  filter_step [label=<<B>6. Filter</B><BR/><FONT POINT-SIZE='8'>(optional)</FONT>>,
               shape=box, fillcolor='#1B4F72', fontcolor='white', style='filled,dashed']

  f_filtered [label=<<B>filtered.fa</B>>,
              shape=note, fillcolor='#D5F5E3', color='#2ECC71', style='filled,dashed']

  // --- MAP ---
  map_step [label=<<B>7. Map</B><BR/><FONT POINT-SIZE='8'>minimap2</FONT>>,
            shape=box, fillcolor='#154360', fontcolor='white']

  f_sam [label=<<B>mapped.sam</B>>,
         shape=note, fillcolor='#D5F5E3', color='#27AE60']

  // --- LIST ---
  list_step [label=<<B>8. List</B><BR/><FONT POINT-SIZE='8'>Count reads</FONT>>,
             shape=box, fillcolor='#154360', fontcolor='white']

  f_detected [label=<<B>detected.list</B>>,
              shape=note, fillcolor='#D5F5E3', color='#27AE60']

  // --- METRICS ---
  metrics [label=<<B>9. Metrics</B><BR/><FONT POINT-SIZE='8'>TP/FP/FN/TN</FONT>>,
           shape=box, fillcolor='#0B5345', fontcolor='white']

  f_results [label=<<B>results.tsv</B>>,
             shape=note, fillcolor='#ABEBC6', color='#1E8449', penwidth=2]
  f_detail [label=<<B>detected_detail.tsv</B>>,
            shape=note, fillcolor='#ABEBC6', color='#1E8449', penwidth=2]
  f_json [label=<<B>results.json</B>>,
          shape=note, fillcolor='#ABEBC6', color='#1E8449', penwidth=2]
  f_coverage [label=<<B>coverage.tsv</B>>,
              shape=note, fillcolor='#ABEBC6', color='#1E8449', penwidth=2]

  // --- REPORT ---
  report [label=<<B>10. Report</B><BR/><FONT POINT-SIZE='8'>R / ggplot2</FONT>>,
          shape=box, fillcolor='#0B5345', fontcolor='white']

  f_report [label=<<B>report.html</B>>,
            shape=note, fillcolor='#82E0AA', color='#1E8449', penwidth=2.5]

  // ==================================================================
  // EDGES: main forward chain (high weight = straight line)
  // ==================================================================

  // Inputs → Prepare
  in_targets -> prepare [weight=3]
  in_distractors -> prepare [weight=3]
  in_sample -> prepare [style=dashed, weight=1]

  // Prepare → its outputs
  prepare -> f_combined [weight=10]
  prepare -> f_weights [weight=5]
  prepare -> f_idlists [weight=3]

  // Prepare outputs → Simulate
  f_combined -> simulate [weight=10]
  f_weights -> simulate [weight=8]

  // Simulate → fragments → Capture
  simulate -> f_fragments [weight=10]
  f_fragments -> capture [weight=10]

  // Probes → Capture (secondary input)
  in_probes -> capture [weight=1]

  // Capture → captured → Enrich → enriched → Sequence
  capture -> f_captured [weight=10]
  f_captured -> enrich [weight=10, style=dashed]
  enrich -> f_enriched [weight=10, style=dashed]
  f_enriched -> sequence [weight=10]

  // Sequence → reads → Filter → filtered → Map
  sequence -> f_reads [weight=10]
  f_reads -> filter_step [weight=10, style=dashed]
  in_host -> filter_step [style=dashed, weight=1]
  filter_step -> f_filtered [weight=10, style=dashed]
  f_filtered -> map_step [weight=10]

  // Reference reuse: combined_reference → Map
  f_combined -> map_step [style=dotted, color='#95A5A6', weight=1,
                          label=<<FONT POINT-SIZE='7' COLOR='#95A5A6'>reference</FONT>>]

  // Map → sam → List → detected → Metrics
  map_step -> f_sam [weight=10]
  f_sam -> list_step [weight=10]
  list_step -> f_detected [weight=10]
  f_detected -> metrics [weight=10]

  // Secondary inputs to Metrics (light dotted)
  f_sam -> metrics [style=dotted, color='#95A5A6', weight=1]
  f_idlists -> metrics [style=dotted, color='#95A5A6', weight=1,
                         label=<<FONT POINT-SIZE='7' COLOR='#95A5A6'>classify</FONT>>]

  // Metrics → outputs
  metrics -> f_results [weight=8]
  metrics -> f_detail [weight=5]
  metrics -> f_json [weight=3]
  metrics -> f_coverage [weight=3]

  // Metrics outputs → Report
  f_results -> report [weight=8]
  f_detail -> report [weight=5]
  f_coverage -> report [style=dashed, weight=2]

  // Report → HTML
  report -> f_report [weight=10]

  // ==================================================================
  // LEGEND
  // ==================================================================
  subgraph cluster_legend {
    label=<<B>Legend</B>>; style=rounded; color='#BDC3C7'; fontsize=10; fontcolor='#95A5A6'
    margin=8
    node [fontsize=9, margin='0.06,0.03']

    leg_step [label='Pipeline step', shape=box, fillcolor='#1A5276', fontcolor='white']
    leg_opt [label='Optional step', shape=box, fillcolor='#21618C', fontcolor='white', style='filled,dashed']
    leg_file [label='Intermediate file', shape=note, fillcolor='#D5F5E3', color='#27AE60']
    leg_input [label='User input', shape=note, fillcolor='#D6EAF8', color='#4A90D9']
    leg_output [label='Final output', shape=note, fillcolor='#ABEBC6', color='#1E8449']

    leg_step -> leg_opt -> leg_file -> leg_input -> leg_output [style=invis]
  }
}
")

save_diagram(pipeline_detailed, "pipeline_detailed.png", width = 4800)
}


# ============================================================================
# Diagram 7: Full BaitBench Tool Workflow — Design → Assess → Simulate → Report
# Three main user-facing workflows shown as phases with explicit data flow.
# ============================================================================
if (7 %in% selected) {
diagram7 <- grViz("
digraph {
  graph [rankdir=LR, fontname='Helvetica', bgcolor='white',
         label=<<B>BaitBench: Complete Tool Workflow</B>>,
         labelloc=t, fontsize=22, pad='0.6,0.4', nodesep=0.45, ranksep=0.9]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box,
        margin='0.15,0.09', penwidth=1.5]
  edge [fontname='Helvetica', fontsize=9, penwidth=1.4]

  // ── USER INPUTS ──────────────────────────────────────────────────────────
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=rounded; color='#BDC3C7'
    fontsize=12; fontcolor='#7F8C8D'; margin=12

    in_targets   [label=<<B>targets.fa</B>>,
                  fillcolor='#D6EAF8', color='#4A90D9', shape=note]
    in_distractors [label=<<B>distractors.fa</B>>,
                    fillcolor='#FADBD8', color='#E74C3C', shape=note]
    in_probes_ext [label=<<B>probes.fa</B>>,
                   fillcolor='#FCF3CF', color='#D4AC0D', shape=note, style='filled,dashed']
    in_sample    [label=<<B>sample.tsv</B>>,
                  fillcolor='#FDEBD0', color='#F39C12', shape=note, style='filled,dashed']

    in_targets -> in_distractors [style=invis]
    in_distractors -> in_probes_ext [style=invis]
    in_probes_ext -> in_sample [style=invis]
  }

  // ── PHASE 1: DESIGN ──────────────────────────────────────────────────────
  subgraph cluster_design {
    label=<<B>Design</B>>
    style=rounded; color='#8E44AD'; fontsize=13; fontcolor='#8E44AD'; margin=15

    build [label=<<B>build-probes</B>>,
           fillcolor='#8E44AD', fontcolor='white']

    probes_fa [label=<<B>probes.fa</B>>,
               fillcolor='#E8DAEF', color='#8E44AD', shape=note]

    build -> probes_fa [color='#8E44AD']
  }

  // ── PHASE 2: ASSESS ──────────────────────────────────────────────────────
  subgraph cluster_assess {
    label=<<B>Assess</B>>
    style=rounded; color='#2980B9'; fontsize=13; fontcolor='#2980B9'; margin=15

    assess [label=<<B>assess-probes</B>>,
            fillcolor='#2980B9', fontcolor='white']

    assess_report [label=<<B>assess_probes.html</B>>,
                   fillcolor='#D6EAF8', color='#2980B9', shape=note]

    assess -> assess_report [color='#2980B9']
  }

  // ── PHASE 3: SIMULATE ────────────────────────────────────────────────────
  subgraph cluster_simulate {
    label=<<B>Simulate</B>>
    style=rounded; color='#1A5276'; fontsize=13; fontcolor='#1A5276'; margin=15

    sim_prepare  [label=<<B>prepare</B>>,     fillcolor='#2C3E50', fontcolor='white']
    sim_simulate [label=<<B>simulate</B>>,    fillcolor='#1A5276', fontcolor='white']
    sim_sequence [label=<<B>sequence</B>>,    fillcolor='#154360', fontcolor='white']
    sim_filter   [label=<<B>filter</B>>,      fillcolor='#154360', fontcolor='white', style='filled,dashed']
    sim_map      [label=<<B>map + list</B>>,  fillcolor='#0D3349', fontcolor='white']
    sim_metrics  [label=<<B>metrics</B>>,     fillcolor='#0B5345', fontcolor='white']

    sim_prepare -> sim_simulate -> sim_sequence -> sim_filter [style=dashed]
    sim_filter -> sim_map [style=dashed]
    sim_sequence -> sim_map [style=dashed]
    sim_map -> sim_metrics
  }

  // ── PHASE 4: REPORT ──────────────────────────────────────────────────────
  subgraph cluster_report {
    label=<<B>Report</B>>
    style=rounded; color='#27AE60'; fontsize=13; fontcolor='#27AE60'; margin=15

    sim_report [label=<<B>report.html</B>>,
                fillcolor='#D5F5E3', color='#27AE60', shape=note, penwidth=2.5]
  }

  // ── EDGES ─────────────────────────────────────────────────────────────────
  // Inputs → Design
  in_targets -> build [color='#8E44AD']

  // Design → Assess
  probes_fa -> assess [color='#2980B9']
  in_targets -> assess [color='#2980B9']

  // External probes → Assess (bypass design)
  in_probes_ext -> assess [style=dashed, color='#2980B9']

  // Inputs → Simulate
  in_targets -> sim_prepare [color='#1A5276']
  in_distractors -> sim_prepare [color='#1A5276']
  in_sample -> sim_prepare [style=dashed, color='#1A5276']

  // Probes → Simulate
  probes_fa -> sim_simulate [color='#1A5276']
  in_probes_ext -> sim_simulate [style=dashed, color='#1A5276']

  // Simulate → Report
  sim_metrics -> sim_report [color='#27AE60', penwidth=2]
}
")

save_diagram(diagram7, "paper_workflow_overview.png", width = 4200)
}


# ============================================================================
# Diagram 8: Thermodynamic Scoring Illustration
# Shows how a probe-reference alignment is converted to a ΔG and Boltzmann score.
# Example: 8-mer probe with one mismatch at position 3.
#   Probe 5'─G─C─A─G─T─C─G─T─3'
#             | | X | | | | |
#   Ref   3'─C─G─C─C─A─G─C─A─5'
#
# WC pairs: positions 1,2,4,5,6,7,8  (n_wc = 7)
# Mismatch: position 3  (A vs C)
# Stacking chains: [1-2] | BREAK | [4-5] [5-6] [6-7] [7-8]
# Initiation: pos 1 = G-C (GC term), pos 8 = T-A (AT term)
# ============================================================================
if (8 %in% selected) {
diagram8 <- grViz("
digraph {
  graph [rankdir=TB, fontname='Helvetica', bgcolor='white',
         label=<<B>Thermodynamic Scoring of a Probe–Reference Alignment</B>>,
         labelloc=t, fontsize=20, pad='0.6,0.4', nodesep=0.6, ranksep=0.65]
  node [fontname='Helvetica', fontsize=11, style=filled, margin='0.15,0.10', penwidth=1.5]
  edge [fontname='Helvetica', fontsize=9, penwidth=1.4, color='#5D6D7E']

  // ── ALIGNMENT TABLE ──────────────────────────────────────────────────────
  alignment [shape=plaintext, label=<
    <TABLE BORDER='0' CELLBORDER='0' CELLSPACING='3' CELLPADDING='5'>

      <TR>
        <TD ALIGN='LEFT'><FONT POINT-SIZE='9' COLOR='#7F8C8D'>Probe  5'</FONT></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>G</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>C</B></TD>
        <TD BGCOLOR='#FADBD8' BORDER='1' WIDTH='32' ALIGN='CENTER'><FONT COLOR='#C0392B'><B>A</B></FONT></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>G</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>T</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>C</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>G</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='32' ALIGN='CENTER'><B>T</B></TD>
        <TD ALIGN='RIGHT'><FONT POINT-SIZE='9' COLOR='#7F8C8D'>3'</FONT></TD>
      </TR>

      <TR>
        <TD></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT COLOR='#C0392B' POINT-SIZE='12'><B>X</B></FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='13'>|</FONT></TD>
        <TD></TD>
      </TR>

      <TR>
        <TD ALIGN='LEFT'><FONT POINT-SIZE='9' COLOR='#7F8C8D'>Ref    3'</FONT></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>G</B></TD>
        <TD BGCOLOR='#FADBD8' BORDER='1' ALIGN='CENTER'><FONT COLOR='#C0392B'><B>C</B></FONT></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>A</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>G</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
        <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>A</B></TD>
        <TD ALIGN='RIGHT'><FONT POINT-SIZE='9' COLOR='#7F8C8D'>5'</FONT></TD>
      </TR>

      <TR>
        <TD></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='8' COLOR='#8E44AD'>init</FONT></TD>
        <TD></TD>
        <TD></TD>
        <TD></TD>
        <TD></TD>
        <TD></TD>
        <TD></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='8' COLOR='#8E44AD'>init</FONT></TD>
        <TD></TD>
      </TR>

      <TR>
        <TD></TD>
        <TD ALIGN='CENTER'
            COLSPAN='2'><FONT POINT-SIZE='8' COLOR='#27AE60'>stack</FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='8' COLOR='#C0392B'>break</FONT></TD>
        <TD ALIGN='CENTER'
            COLSPAN='5'><FONT POINT-SIZE='8' COLOR='#27AE60'>stack  stack  stack  stack</FONT></TD>
        <TD></TD>
      </TR>

    </TABLE>
  >]

  // ── KEY ──────────────────────────────────────────────────────────────────
  legend_align [shape=plaintext, label=<
    <TABLE BORDER='1' CELLBORDER='0' CELLSPACING='0' CELLPADDING='5' COLOR='#BDC3C7'>
      <TR><TD BGCOLOR='#D5F5E3' WIDTH='18'> </TD>
          <TD ALIGN='LEFT'><FONT POINT-SIZE='9'>Watson-Crick pair</FONT></TD></TR>
      <TR><TD BGCOLOR='#FADBD8' WIDTH='18'> </TD>
          <TD ALIGN='LEFT'><FONT POINT-SIZE='9'>Mismatch (breaks stacking chain)</FONT></TD></TR>
    </TABLE>
  >]

  // ── THREE CONTRIBUTIONS ───────────────────────────────────────────────────
  stacking [shape=box, fillcolor='#D5F5E3', color='#27AE60',
            label=<<B>1. NN Stacking</B><BR/>
            <FONT POINT-SIZE='9'>For each consecutive WC step:<BR/>
            add ΔH and ΔS from SantaLucia (1998)<BR/>
            Table 2 (10 dinucleotide parameters).<BR/><BR/>
            Mismatch at pos 3 breaks chain;<BR/>
            5 stacking steps contribute here.</FONT>>]

  initiation [shape=box, fillcolor='#E8DAEF', color='#8E44AD',
              label=<<B>2. Initiation Terms</B><BR/>
              <FONT POINT-SIZE='9'>First WC pair (pos 1, G-C):<BR/>
              +0.1 kcal/mol ΔH,  -2.8 cal/mol/K ΔS<BR/><BR/>
              Last WC pair (pos 8, T-A):<BR/>
              +2.3 kcal/mol ΔH,  +4.1 cal/mol/K ΔS<BR/><BR/>
              (SantaLucia 1998 Table 2)</FONT>>]

  salt [shape=box, fillcolor='#FEF9E7', color='#D4AC0D',
        label=<<B>3. Salt Correction</B><BR/>
        <FONT POINT-SIZE='9'>n_wc = 7 Watson-Crick pairs<BR/><BR/>
        ΔS += 0.368 × (n_wc − 1) × ln([Na<SUP>+</SUP>])<BR/>
        e.g. at 50 mM:  ln(0.05) ≈ −3.0<BR/>
        correction ≈ −6.6 cal/mol/K per probe<BR/><BR/>
        (Owczarzy et al. 1997)<BR/>
        At 1 M Na<SUP>+</SUP>: correction = 0</FONT>>]

  // ── ΔG CALCULATION ───────────────────────────────────────────────────────
  delta_g_box [shape=box, fillcolor='#D6EAF8', color='#2980B9', penwidth=2,
               label=<<B>ΔG (kcal/mol)</B><BR/>
               <FONT POINT-SIZE='10'>ΔG = ΔH<SUB>stack+init</SUB> − T × (ΔS<SUB>stack+init+salt</SUB> / 1000)</FONT><BR/>
               <FONT POINT-SIZE='9'>T in Kelvin  |  more negative = more stable</FONT>>]

  // ── BOLTZMANN SCORE ───────────────────────────────────────────────────────
  boltzmann [shape=box, fillcolor='#1A5276', fontcolor='white', penwidth=2,
             label=<<B>Boltzmann Binding Score</B><BR/>
             <FONT POINT-SIZE='10'>score = exp( −ΔG / (R × T) )</FONT><BR/>
             <FONT POINT-SIZE='9'>R = 1.987 × 10<SUP>−3</SUP> kcal/(mol·K)<BR/>
             Higher score → stronger predicted binding</FONT>>]

  // ── SAMPLING WEIGHT ───────────────────────────────────────────────────────
  weight [shape=box, fillcolor='#0B5345', fontcolor='white',
          label=<<B>Sampling Weight</B><BR/>
          <FONT POINT-SIZE='9'>weight = Boltzmann score × sequence weight<BR/>
          Used in two-level multinomial sampling<BR/>
          to bias fragments toward high-affinity sites</FONT>>]

  // ── EDGES ─────────────────────────────────────────────────────────────────
  alignment -> stacking   [label=<<FONT POINT-SIZE='8'>WC pair stream</FONT>>]
  alignment -> initiation [label=<<FONT POINT-SIZE='8'>terminal WC pairs</FONT>>]
  alignment -> salt       [label=<<FONT POINT-SIZE='8'>n_wc count</FONT>>]
  alignment -> legend_align [style=invis]

  stacking   -> delta_g_box [label=<<FONT POINT-SIZE='8'>sum ΔH, ΔS</FONT>>]
  initiation -> delta_g_box [label=<<FONT POINT-SIZE='8'>add ΔH, ΔS</FONT>>]
  salt       -> delta_g_box [label=<<FONT POINT-SIZE='8'>correct ΔS</FONT>>]

  delta_g_box -> boltzmann [penwidth=2]
  boltzmann   -> weight    [penwidth=2]
}
")

save_diagram(diagram8, "paper_thermodynamic_scoring.png", width = 3600)
}


# ============================================================================
# Diagram 9: How a Fragment Gets Selected
# Shows the two arms by which a fragment reaches fragments.fa.
# ============================================================================
# Layout notes: <I> around a single subscript character makes graphviz
# mis-measure the run, producing huge gaps and breaking line centring. Use
# plain _subscript notation instead. Likewise never mix inline <B> with plain
# text on the same line.
if (9 %in% selected) {
diagram9 <- grViz("
digraph {
  graph [rankdir=TB, fontname='Helvetica', bgcolor='white',
         label=<<B>How a Fragment Gets Selected</B><BR/>
         <FONT POINT-SIZE='14'>Two arms: probe-mediated capture vs non-specific bleed-through</FONT>>,
         labelloc=t, fontsize=22, pad='0.6,0.4', nodesep=0.55, ranksep=0.7]
  node [fontname='Helvetica', fontsize=11, style=filled, margin='0.15,0.10', penwidth=1.5]
  edge [fontname='Helvetica', fontsize=9, penwidth=1.4, color='#5D6D7E']

  // ── INPUTS ───────────────────────────────────────────────────────────────
  inputs [shape=plaintext, label=<
    <TABLE BORDER='1' CELLBORDER='0' CELLSPACING='0' CELLPADDING='7'
           COLOR='#BDC3C7' BGCOLOR='#F2F3F4'>
      <TR>
        <TD COLSPAN='3' ALIGN='CENTER' BGCOLOR='#D6EAF8'>
          <B>Simulation inputs</B>
        </TD>
      </TR>
      <TR>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>N total fragments<BR/><B>--num-fragments</B></FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>capture fraction f<BR/><B>--capture-fraction</B></FONT></TD>
        <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>probe alignments<BR/>(SAM, TNN-scored)</FONT></TD>
      </TR>
    </TABLE>
  >]

  // ── SPLIT ────────────────────────────────────────────────────────────────
  split [shape=diamond, fillcolor='#34495E', fontcolor='white', penwidth=2,
         label=<<B>Split N fragments</B>>]

  n_capture   [shape=box, fillcolor='#1A5276', fontcolor='white',
               label=<<B>N × f fragments</B><BR/><FONT POINT-SIZE='9'>probe-biased (captured pool)</FONT>>]
  n_background [shape=box, fillcolor='#6D7C8A', fontcolor='white',
                label=<<B>N × (1 − f) fragments</B><BR/><FONT POINT-SIZE='9'>background (uncaptured)</FONT>>]

  // ── ARM HEADINGS (above each step box) ───────────────────────────────────
  cap_title [shape=plaintext, fillcolor='white',
             label=<<FONT POINT-SIZE='15' COLOR='#1A5276'><B>Probe-Biased Sampling</B></FONT><BR/>
             <FONT POINT-SIZE='11' COLOR='#1A5276'>Two-level weighted multinomial</FONT>>]

  bg_title [shape=plaintext, fillcolor='white',
            label=<<FONT POINT-SIZE='15' COLOR='#5D6D7E'><B>Background Sampling</B></FONT><BR/>
            <FONT POINT-SIZE='11' COLOR='#5D6D7E'>Uniform by sequence weight</FONT>>]

  // ── PROBE-BIASED BRANCH (left) ────────────────────────────────────────────
  subgraph cluster_capture {
    label=\"\"; style=rounded; color='#1A5276'; margin=14; penwidth=2

    step1 [shape=box, fillcolor='#2980B9', fontcolor='white',
           label=<<B>Select probe</B><BR/>
           <FONT POINT-SIZE='9'>Uniform over probes<BR/>
           with at least 1 alignment hit<BR/><BR/>
           P(probe_i) = 1 / number of probes with hits</FONT>>]

    step2 [shape=box, fillcolor='#1A5276', fontcolor='white',
           label=<<B>Select alignment hit</B><BR/>
           <FONT POINT-SIZE='9'>Weighted by Boltzmann score × sequence weight<BR/><BR/>
           P(hit_j) proportional to exp(−ΔG_j / RT) × w_seq<BR/><BR/>
           High-affinity sites preferentially sampled</FONT>>]

    step3 [shape=box, fillcolor='#154360', fontcolor='white',
           label=<<B>Fragment center</B><BR/>
           <FONT POINT-SIZE='9'>center = hit_center + Uniform(−L/4, +L/4)<BR/>
           where L = fragment length mean<BR/><BR/>
           Jitter models probe offset variability</FONT>>]

    step4 [shape=box, fillcolor='#0D3349', fontcolor='white',
           label=<<B>Fragment length</B><BR/>
           <FONT POINT-SIZE='9'>length ~ TruncNormal(mean, sd)<BR/>
           clamped to [min, max]<BR/><BR/>
           Extract sequence from reference</FONT>>]

    step1 -> step2 -> step3 -> step4 [color='#2980B9', penwidth=1.8]
  }

  // ── BACKGROUND BRANCH (right) ─────────────────────────────────────────────
  subgraph cluster_background {
    label=\"\"; style=rounded; color='#6D7C8A'; margin=14; penwidth=2

    bg1 [shape=box, fillcolor='#6D7C8A', fontcolor='white',
         label=<<B>Select sequence</B><BR/>
         <FONT POINT-SIZE='9'>P(seq_i) proportional to weight_i × length_i<BR/><BR/>
         Sequence weights from weights.txt<BR/>
         (sample targets, distractors)</FONT>>]

    bg2 [shape=box, fillcolor='#5D6D7E', fontcolor='white',
         label=<<B>Generate fragment</B><BR/>
         <FONT POINT-SIZE='9'>Position ~ Uniform(0, seq_len − frag_len)<BR/>
         length ~ TruncNormal(mean, sd)<BR/><BR/>
         No probe-site bias</FONT>>]

    bg1 -> bg2 [color='#6D7C8A', penwidth=1.8]
  }

  // ── OUTPUT ───────────────────────────────────────────────────────────────
  merge [shape=diamond, fillcolor='#0B5345', fontcolor='white', penwidth=2,
         label=<<B>Combine both pools</B>>]

  output [shape=folder, fillcolor='#D5F5E3', color='#27AE60', penwidth=2.5,
          label=<<B>fragments.fa</B><BR/>
          <FONT POINT-SIZE='9'>N fragments total<BR/>
          Names encode source: {seq_id}_fragment_{n}<BR/>
          start=X length=Y</FONT>>]

  // ── ANNOTATION ───────────────────────────────────────────────────────────
  note_cf [shape=note, fillcolor='#FEF9E7', color='#D4AC0D',
           label=<<FONT POINT-SIZE='9'><I>capture-fraction models overall<BR/>
           pull-down efficiency (wash stringency,<BR/>
           hybridization time). Per-site differential<BR/>
           enrichment is handled by Boltzmann<BR/>
           weighting at the hit-selection step.</I></FONT>>]

  // ── EDGES ─────────────────────────────────────────────────────────────────
  inputs -> split [penwidth=2]
  split -> n_capture [color='#1A5276', penwidth=2, minlen=2,
                      label=<<FONT POINT-SIZE='26' COLOR='#1A5276'><B>&nbsp;&nbsp;f&nbsp;&nbsp;</B></FONT>>]
  split -> n_background [color='#6D7C8A', penwidth=2, minlen=2,
                         label=<<FONT POINT-SIZE='26' COLOR='#5D6D7E'><B>&nbsp;&nbsp;1 − f&nbsp;&nbsp;</B></FONT>>]

  n_capture -> cap_title [color='#1A5276', penwidth=1.8, arrowhead=none]
  cap_title -> step1 [color='#1A5276', penwidth=1.8]

  n_background -> bg_title [color='#6D7C8A', penwidth=1.8, arrowhead=none]
  bg_title -> bg1 [color='#6D7C8A', penwidth=1.8]

  step4 -> merge [color='#1A5276', penwidth=1.8]
  bg2   -> merge [color='#6D7C8A', penwidth=1.8]

  merge -> output [penwidth=2.5, color='#27AE60']

  note_cf -> split [style=invis]
}
")

save_diagram(diagram9, "paper_fragment_sampling.png", width = 3800)
}



# ============================================================================
# Diagram 10: build-probes pipeline steps
# ============================================================================
if (10 %in% selected) {
diagram10 <- grViz("
digraph {
  graph [rankdir=TB, fontname='Helvetica', bgcolor='white',
         label=<<B>baitbench build-probes</B>>,
         labelloc=t, fontsize=20, pad=0.6, nodesep=0.55, ranksep=0.7]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box,
        margin='0.18,0.10']
  edge [fontname='Helvetica', fontsize=10]

  // ── INPUT ─────────────────────────────────────────────────────────────────
  input [label=<<B>targets.fa</B>>,
         fillcolor='#D6EAF8', color='#4A90D9', penwidth=2, shape=folder]

  // ── STEPS ─────────────────────────────────────────────────────────────────
  s1 [label=<<B>N-content filter</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  s2 [label=<<B>Collapse redundant targets</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  s3 [label=<<B>Length filter</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  // ── method branch ────────────────────────────────────────────────────────
  s4_choice [label=<<B>Build probes</B>>,
             fillcolor='#D5F5E3', color='#1E8449', penwidth=2, shape=diamond,
             margin='0.25,0.12']

  s4_tile   [label=<<B>tile</B>>,        fillcolor='#EAFAF1', color='#1E8449', penwidth=1.5]
  s4_catch  [label=<<B>catch-lite</B>>,  fillcolor='#EAFAF1', color='#1E8449', penwidth=1.5]
  s4_syotti [label=<<B>syotti-lite</B>>, fillcolor='#EAFAF1', color='#1E8449', penwidth=1.5]
  s4_ext    [label=<<B>catch</B>>,       fillcolor='#EAFAF1', color='#1E8449', penwidth=1.5]

  s4_merge [label='', shape=point, width=0.15, fillcolor='#1E8449', color='#1E8449']

  s5 [label=<<B>GC-content filter</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  s6 [label=<<B>Complexity filter</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  s7 [label=<<B>Deduplicate</B>>,
      fillcolor='#FDEBD0', color='#CA6F1E', penwidth=1.5]

  // ── OUTPUT ────────────────────────────────────────────────────────────────
  output [label=<<B>probes.fa</B>>,
          fillcolor='#D5F5E3', color='#1E8449', penwidth=2, shape=folder]

  // ── EDGES ─────────────────────────────────────────────────────────────────
  input -> s1 -> s2 -> s3 -> s4_choice [penwidth=2]

  s4_choice -> s4_tile   [label='tile',       color='#1E8449']
  s4_choice -> s4_catch  [label='catch-lite', color='#1E8449']
  s4_choice -> s4_syotti [label='syotti-lite',color='#1E8449']
  s4_choice -> s4_ext    [label='catch',      color='#1E8449']

  s4_tile   -> s4_merge [color='#1E8449']
  s4_catch  -> s4_merge [color='#1E8449']
  s4_syotti -> s4_merge [color='#1E8449']
  s4_ext    -> s4_merge [color='#1E8449']

  s4_merge -> s5 -> s6 -> s7 -> output [penwidth=2]
}
")

save_diagram(diagram10, "paper_build_probes.png", width = 3600)
}


# ============================================================================
# Diagram 11: assess-probes pipeline steps
# ============================================================================
if (11 %in% selected) {
diagram11 <- grViz("
digraph {
  graph [rankdir=TB, fontname='Helvetica', bgcolor='white',
         label=<<B>baitbench assess-probes</B>>,
         labelloc=t, fontsize=20, pad=0.6, nodesep=0.55, ranksep=0.7]
  node [fontname='Helvetica', fontsize=12, style=filled, shape=box,
        margin='0.18,0.10']
  edge [fontname='Helvetica', fontsize=10]

  // ── INPUTS ────────────────────────────────────────────────────────────────
  subgraph cluster_inputs {
    label=<<B>Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=13;
    fontcolor='#7F8C8D'

    targets [label=<<B>targets.fa</B>>,
             fillcolor='#D6EAF8', color='#4A90D9', penwidth=2, shape=folder]

    probes [label=<<B>probes.fa</B>>,
            fillcolor='#D5F5E3', color='#1E8449', penwidth=2, shape=folder]

    genomes_opt [label=<<B>genomes.fa</B>>,
                 fillcolor='#F2F3F4', color='#95A5A6', penwidth=1, shape=folder]

    build_stats [label=<<B>build_stats.tsv</B>>,
                 fillcolor='#F2F3F4', color='#95A5A6', penwidth=1, shape=folder]
  }

  // ── COVERAGE ──────────────────────────────────────────────────────────────
  s1a [label=<<B>Probe coverage analysis</B>>,
       fillcolor='#EBF5FB', color='#2E86C1', penwidth=2]

  s1b [label=<<B>Individual target coverage</B>>,
       fillcolor='#EBF5FB', color='#2E86C1', penwidth=2]

  s1c [label=<<B>Gap details</B>>,
       fillcolor='#EBF5FB', color='#2E86C1', penwidth=2]

  // ── CROSS-REACTIVITY ──────────────────────────────────────────────────────
  s2a [label=<<B>Self cross-reactivity</B>>,
       fillcolor='#FDEDEC', color='#C0392B', penwidth=2]

  s2b [label=<<B>Genome cross-reactivity</B>>,
       fillcolor='#FDEDEC', color='#C0392B', penwidth=1.5, style='filled,dashed']

  // ── REPORT ────────────────────────────────────────────────────────────────
  report [label=<<B>Combined HTML Report</B>>,
          fillcolor='#F9F0FF', color='#7D3C98', penwidth=2, shape=box]

  output [label=<<B>assess_probes_report.html</B>>,
          fillcolor='#F3E5F5', color='#7D3C98', penwidth=2, shape=folder]

  // ── EDGES ─────────────────────────────────────────────────────────────────
  targets -> s1a [penwidth=2]
  probes  -> s1a [penwidth=2]
  probes  -> s2a [penwidth=2]

  targets -> s1b [penwidth=2]
  probes  -> s1b [penwidth=2]
  s1a -> s1c [color='#2E86C1', penwidth=2]
  s1b -> s1c [color='#2E86C1', penwidth=2]

  genomes_opt -> s2b [style=dashed, color='#95A5A6', penwidth=1.2]
  probes -> s2b [style=dashed, color='#C0392B', penwidth=1.2]

  build_stats -> report [style=dashed, color='#95A5A6', penwidth=1.2]

  s1a -> report [penwidth=2, color='#2E86C1']
  s1b -> report [penwidth=2, color='#2E86C1']
  s1c -> report [penwidth=2, color='#2E86C1']
  s2a -> report [penwidth=2, color='#C0392B']
  s2b -> report [penwidth=1.2, style=dashed, color='#C0392B']

  report -> output [penwidth=2.5, color='#7D3C98']
}
")

save_diagram(diagram11, "paper_assess_probes.png", width = 3600)
}



# ============================================================================
# Diagram 12: Full Simulation Pipeline — Genome Mode, Full Detail
# Deliberately dense: merges mode 4 (--genomes + --sample) prepare, probe
# alignment, thermodynamic scoring (diagram 8), fragment sampling (diagram 9),
# every read-simulator option, filter, map/list, genome-aware metrics, and
# report generation into one diagram. Not meant to be read at a glance.
# ============================================================================
if (12 %in% selected) {
diagram12 <- grViz("
digraph {
  graph [rankdir=TB, fontname='Helvetica', bgcolor='white',
         label=<<B>BaitBench: Full Simulation Pipeline — Genome Mode</B><BR/>
         <FONT POINT-SIZE='13'>Mode 4 (--genomes + --sample) &middot; thermodynamic capture model &middot; every stage in detail</FONT>>,
         labelloc=t, fontsize=24, pad='0.7,0.5', nodesep=0.45, ranksep=0.75]
  node [fontname='Helvetica', fontsize=11, style=filled, shape=box, margin='0.14,0.08']
  edge [fontname='Helvetica', fontsize=9]

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 1 — USER INPUTS
  // ══════════════════════════════════════════════════════════════════════
  subgraph cluster_inputs {
    label=<<B>User Inputs</B>>; style=dashed; color='#7F8C8D'; fontsize=14; fontcolor='#7F8C8D'

    in_genomes [label=<<B>genomes.fa</B><BR/><FONT POINT-SIZE='10'>Full genome sequences<BR/>IDs: G1, G2, G3</FONT>>,
                fillcolor='#E8DAEF', color='#8E44AD', penwidth=2]
    in_targets [label=<<B>targets.fa</B><BR/><FONT POINT-SIZE='10'>Probe target subsequences<BR/>IDs: G1|16S, G2|ompB, G3|gltA</FONT>>,
                fillcolor='#D6EAF8', color='#4A90D9', penwidth=2]
    in_distractors [label=<<B>distractors.fa</B><BR/><FONT POINT-SIZE='10'>Background sequences<BR/>IDs: D1, D2, ... Dn</FONT>>,
                    fillcolor='#FADBD8', color='#E74C3C', penwidth=2]
    in_probes [label=<<B>probes.fa</B><BR/><FONT POINT-SIZE='10'>Capture probe sequences<BR/>(from build-probes or user-supplied)</FONT>>,
               fillcolor='#FCF3CF', color='#D4AC0D', penwidth=2]
    in_sample [label=<<B>--sample manifest.tsv</B><BR/><FONT POINT-SIZE='10'>G1  1.0<BR/>G2  3.0<BR/>(G3 is NOT in sample)</FONT>>,
               fillcolor='#FDEBD0', color='#F39C12', penwidth=2, shape=note]
    in_stmap [label=<<B>--sample-target-map</B><BR/><FONT POINT-SIZE='10'>(optional explicit TSV)<BR/>G1  G1|16S<BR/>G2  G2|ompB<BR/>G3  G3|gltA<BR/><I>Or: auto-linked by prefix</I></FONT>>,
              fillcolor='#FDEBD0', color='#E67E22', penwidth=2, shape=note]
    in_params [label=<<B>Global Parameters</B><BR/>
               <FONT POINT-SIZE='9'>--distractor-fraction 0.9 (or --ct)<BR/>
               --hybridization-temperature 65<BR/>
               --salt-concentration 50<BR/>
               --fragment-length-mean/min/max<BR/>
               --capture-fraction<BR/>
               --num-fragments<BR/>
               --seed</FONT>>,
               fillcolor='#F2F3F4', color='#95A5A6', shape=note, penwidth=1]
  }

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 2 — PREPARE
  // ══════════════════════════════════════════════════════════════════════
  prep_step [label=<<B>baitbench prepare</B><BR/><FONT POINT-SIZE='10'>Validate sample IDs in genomes<BR/>Resolve genome &rarr; target maps<BR/>Build two references<BR/>Generate weights</FONT>>,
             fillcolor='#34495E', fontcolor=white, shape=box, penwidth=2]

  subgraph cluster_prepare_out {
    label=<<B>Prepare Outputs</B>>; style=dashed; color='#27AE60'; fontsize=14; fontcolor='#27AE60'

    prep_combined [label=<<B>combined_reference.fa</B><BR/><FONT POINT-SIZE='10'>genomes.fa + distractors.fa<BR/>For fragment generation</FONT>>,
                   fillcolor='#D5F5E3', color='#27AE60', penwidth=2]
    prep_mapping [label=<<B>mapping_reference.fa</B><BR/><FONT POINT-SIZE='10'>targets.fa + distractors.fa<BR/>For read mapping (later)</FONT>>,
                  fillcolor='#D5F5E3', color='#27AE60', penwidth=2]
    prep_weights [label=<<B>weights.txt</B><BR/><FONT POINT-SIZE='10'>G1  1.0   (sample)<BR/>G2  3.0   (sample, high weight)<BR/><FONT COLOR='#95A5A6'>G3  0.0   (non-sample, no frags!)</FONT><BR/>D1  0.036 (distractor weight)<BR/>...</FONT>>,
                  fillcolor='#D5F5E3', color='#2ECC71', penwidth=2]
    prep_stmap [label=<<B>sample_target_map.txt</B><BR/><FONT POINT-SIZE='10'>G1 &rarr; G1|16S<BR/>G2 &rarr; G2|ompB<BR/>G3 &rarr; G3|gltA<BR/>(used by metrics step)</FONT>>,
                fillcolor='#FDEBD0', color='#E67E22', penwidth=2]
    prep_idlists [label=<<B>ID Lists</B><BR/><FONT POINT-SIZE='10'>genomes.txt, targets.txt<BR/>distractors.txt<BR/>sample.txt (subset!)</FONT>>,
                  fillcolor='#D1F2EB', color='#1ABC9C', penwidth=2]
  }

  in_genomes -> prep_step
  in_targets -> prep_step
  in_distractors -> prep_step
  in_sample -> prep_step
  in_stmap -> prep_step [style=dashed, label=<<FONT POINT-SIZE='9'>optional</FONT>>]
  in_params -> prep_step [style=dashed, color='#95A5A6', arrowhead=none]

  prep_step -> prep_combined
  prep_step -> prep_mapping
  prep_step -> prep_weights
  prep_step -> prep_stmap
  prep_step -> prep_idlists

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 3 — SIMULATE: probe alignment entry
  // ══════════════════════════════════════════════════════════════════════
  sim_align [label=<<B>Probe Alignment (rammap)</B><BR/>
             <FONT POINT-SIZE='10'>probes.fa vs combined_reference.fa<BR/>
             &rarr; fragments.probe_hits.sam</FONT>>,
             fillcolor='#2C3E50', fontcolor=white, shape=box, penwidth=2]

  in_probes -> sim_align
  prep_combined -> sim_align

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 4 — THERMODYNAMIC SCORING  (SantaLucia 1998 NN model)
  // ══════════════════════════════════════════════════════════════════════
  subgraph cluster_thermo {
    label=<<B>Thermodynamic Scoring — per probe/reference alignment</B>>;
    style=rounded; color='#2980B9'; fontsize=13; fontcolor='#2980B9'; margin=16

    td_alignment [shape=plaintext, label=<
      <TABLE BORDER='0' CELLBORDER='0' CELLSPACING='2' CELLPADDING='4'>
        <TR>
          <TD ALIGN='LEFT'><FONT POINT-SIZE='8' COLOR='#7F8C8D'>Probe 5'</FONT></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>G</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>C</B></TD>
          <TD BGCOLOR='#FADBD8' BORDER='1' WIDTH='24' ALIGN='CENTER'><FONT COLOR='#C0392B'><B>A</B></FONT></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>G</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>T</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>C</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>G</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' WIDTH='24' ALIGN='CENTER'><B>T</B></TD>
          <TD ALIGN='RIGHT'><FONT POINT-SIZE='8' COLOR='#7F8C8D'>3'</FONT></TD>
        </TR>
        <TR>
          <TD></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT COLOR='#C0392B' POINT-SIZE='10'><B>X</B></FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='11'>|</FONT></TD>
          <TD></TD>
        </TR>
        <TR>
          <TD ALIGN='LEFT'><FONT POINT-SIZE='8' COLOR='#7F8C8D'>Ref   3'</FONT></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>G</B></TD>
          <TD BGCOLOR='#FADBD8' BORDER='1' ALIGN='CENTER'><FONT COLOR='#C0392B'><B>C</B></FONT></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>A</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>G</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>C</B></TD>
          <TD BGCOLOR='#D5F5E3' BORDER='1' ALIGN='CENTER'><B>A</B></TD>
          <TD ALIGN='RIGHT'><FONT POINT-SIZE='8' COLOR='#7F8C8D'>5'</FONT></TD>
        </TR>
      </TABLE>
    >]

    td_stacking [shape=box, fillcolor='#D5F5E3', color='#27AE60',
                 label=<<B>1. NN Stacking</B><BR/>
                 <FONT POINT-SIZE='9'>Sum &#916;H, &#916;S per consecutive<BR/>
                 WC step (SantaLucia 1998, Table 2:<BR/>
                 10 dinucleotide parameters).<BR/>
                 Mismatch breaks the chain.</FONT>>]

    td_initiation [shape=box, fillcolor='#E8DAEF', color='#8E44AD',
                   label=<<B>2. Initiation Terms</B><BR/>
                   <FONT POINT-SIZE='9'>First + last WC pair each add<BR/>
                   a fixed &#916;H/&#916;S initiation term<BR/>
                   (SantaLucia 1998, Table 2)</FONT>>]

    td_salt [shape=box, fillcolor='#FEF9E7', color='#D4AC0D',
             label=<<B>3. Salt Correction</B><BR/>
             <FONT POINT-SIZE='9'>&#916;S += 0.368 &#215; (n_wc &#8722; 1) &#215; ln([Na+])<BR/>
             (Owczarzy et al. 1997)<BR/>
             n_wc = count of WC pairs<BR/>
             At 1 M Na+: correction = 0</FONT>>]

    td_deltag [shape=box, fillcolor='#D6EAF8', color='#2980B9', penwidth=2,
               label=<<B>&#916;G (kcal/mol)</B><BR/>
               <FONT POINT-SIZE='9'>&#916;G = &#916;H &#8722; T &#215; (&#916;S / 1000)<BR/>
               more negative = more stable</FONT>>]

    td_boltzmann [shape=box, fillcolor='#1A5276', fontcolor='white', penwidth=2,
                  label=<<B>Boltzmann Binding Score</B><BR/>
                  <FONT POINT-SIZE='9'>score = exp(&#8722;&#916;G / (R &#215; T))<BR/>
                  R = 1.987&#215;10<SUP>&#8722;3</SUP> kcal/(mol&middot;K)</FONT>>]

    td_weight [shape=box, fillcolor='#0B5345', fontcolor='white', penwidth=2,
               label=<<B>Per-Hit Sampling Weight</B><BR/>
               <FONT POINT-SIZE='9'>weight = Boltzmann score &#215; sequence weight<BR/>
               feeds probe-biased sampling below</FONT>>]

    td_alignment -> td_stacking   [label=<<FONT POINT-SIZE='8'>WC pairs</FONT>>]
    td_alignment -> td_initiation [label=<<FONT POINT-SIZE='8'>terminal pairs</FONT>>]
    td_alignment -> td_salt       [label=<<FONT POINT-SIZE='8'>n_wc</FONT>>]
    td_stacking   -> td_deltag [label=<<FONT POINT-SIZE='8'>&#916;H, &#916;S</FONT>>]
    td_initiation -> td_deltag [label=<<FONT POINT-SIZE='8'>&#916;H, &#916;S</FONT>>]
    td_salt       -> td_deltag [label=<<FONT POINT-SIZE='8'>&#916;S corr.</FONT>>]
    td_deltag -> td_boltzmann [penwidth=1.8]
    td_boltzmann -> td_weight [penwidth=1.8]
  }

  sim_align -> td_alignment [penwidth=2, label=<<FONT POINT-SIZE='8'>per hit</FONT>>]
  sim_align -> fs_inputs [penwidth=2, style=dashed, color='#5D6D7E',
                          label=<<FONT POINT-SIZE='8'>--num-fragments, --capture-fraction</FONT>>]
  prep_weights -> td_weight [style=dashed, color='#95A5A6']

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 5 — FRAGMENT SAMPLING  (two-level multinomial)
  // ══════════════════════════════════════════════════════════════════════
  subgraph cluster_fragsamp {
    label=<<B>Fragment Sampling — two arms per fragment</B>>;
    style=rounded; color='#1A5276'; fontsize=13; fontcolor='#1A5276'; margin=16

    fs_inputs [shape=plaintext, label=<
      <TABLE BORDER='1' CELLBORDER='0' CELLSPACING='0' CELLPADDING='6' COLOR='#BDC3C7' BGCOLOR='#F2F3F4'>
        <TR><TD COLSPAN='3' ALIGN='CENTER' BGCOLOR='#D6EAF8'><B>Simulation inputs</B></TD></TR>
        <TR>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>N total fragments<BR/><B>--num-fragments</B></FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>capture fraction f<BR/><B>--capture-fraction</B></FONT></TD>
          <TD ALIGN='CENTER'><FONT POINT-SIZE='9'>weighted probe hits<BR/>(from thermo scoring)</FONT></TD>
        </TR>
      </TABLE>
    >]

    fs_split [shape=diamond, fillcolor='#34495E', fontcolor='white', penwidth=2,
              label=<<B>Split N fragments</B>>]

    fs_ncapture [shape=box, fillcolor='#1A5276', fontcolor='white',
                 label=<<B>N &#215; f fragments</B><BR/><FONT POINT-SIZE='9'>probe-biased (captured pool)</FONT>>]
    fs_nbackground [shape=box, fillcolor='#6D7C8A', fontcolor='white',
                    label=<<B>N &#215; (1 &#8722; f) fragments</B><BR/><FONT POINT-SIZE='9'>background (uncaptured)</FONT>>]

    fs_step1 [shape=box, fillcolor='#2980B9', fontcolor='white',
              label=<<B>Select probe</B><BR/><FONT POINT-SIZE='9'>Uniform over probes with<BR/>at least 1 alignment hit</FONT>>]
    fs_step2 [shape=box, fillcolor='#1A5276', fontcolor='white',
              label=<<B>Select alignment hit</B><BR/><FONT POINT-SIZE='9'>P(hit) proportional to<BR/>exp(&#8722;&#916;G/RT) &#215; w_seq<BR/>high-affinity sites favored</FONT>>]
    fs_step3 [shape=box, fillcolor='#154360', fontcolor='white',
              label=<<B>Fragment center</B><BR/><FONT POINT-SIZE='9'>center = hit_center<BR/>+ Uniform(&#8722;L/4, +L/4)</FONT>>]
    fs_step4 [shape=box, fillcolor='#0D3349', fontcolor='white',
              label=<<B>Fragment length</B><BR/><FONT POINT-SIZE='9'>length ~ TruncNormal(mean, sd)<BR/>clamped to [min, max]</FONT>>]

    fs_bg1 [shape=box, fillcolor='#6D7C8A', fontcolor='white',
            label=<<B>Select sequence</B><BR/><FONT POINT-SIZE='9'>P(seq) proportional to<BR/>weight &#215; length (weights.txt)</FONT>>]
    fs_bg2 [shape=box, fillcolor='#5D6D7E', fontcolor='white',
            label=<<B>Generate fragment</B><BR/><FONT POINT-SIZE='9'>Uniform position, TruncNormal length<BR/>no probe-site bias</FONT>>]

    fs_merge [shape=diamond, fillcolor='#0B5345', fontcolor='white', penwidth=2,
              label=<<B>Combine both pools</B>>]

    fs_inputs -> fs_split [penwidth=2]
    fs_split -> fs_ncapture    [color='#1A5276', penwidth=2, label=<<FONT POINT-SIZE='16' COLOR='#1A5276'><B>f</B></FONT>>]
    fs_split -> fs_nbackground [color='#6D7C8A', penwidth=2, label=<<FONT POINT-SIZE='16' COLOR='#5D6D7E'><B>1&#8722;f</B></FONT>>]
    fs_ncapture -> fs_step1 [color='#1A5276']
    fs_step1 -> fs_step2 -> fs_step3 -> fs_step4 [color='#2980B9']
    fs_nbackground -> fs_bg1 [color='#6D7C8A']
    fs_bg1 -> fs_bg2 [color='#6D7C8A']
    fs_step4 -> fs_merge [color='#1A5276']
    fs_bg2 -> fs_merge [color='#6D7C8A']
  }

  td_weight -> fs_inputs [penwidth=2.5, color='#2980B9', constraint=false,
                          label=<<FONT POINT-SIZE='9' COLOR='#2980B9'>per-hit weight</FONT>>]
  in_params -> fs_inputs [style=dashed, color='#95A5A6', arrowhead=none]
  prep_weights -> fs_bg1 [style=dashed, color='#95A5A6']

  fs_output [label=<<B>fragments.fa</B><BR/><FONT POINT-SIZE='10'>N fragments total<BR/>Names encode source:<BR/>{seq_id}_fragment_{n}</FONT>>,
             fillcolor='#D5F5E3', color='#27AE60', penwidth=2.5, shape=folder]

  fs_merge -> fs_output [penwidth=2.5, color='#27AE60']

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 6 — SEQUENCE: read simulator choice
  // ══════════════════════════════════════════════════════════════════════
  subgraph cluster_sequence {
    label=<<B>baitbench sequence — read simulator options</B>>;
    style=rounded; color='#117864'; fontsize=13; fontcolor='#117864'; margin=16

    seq_sample [label=<<B>Optional subsample</B><BR/><FONT POINT-SIZE='9'>--num-sequences N<BR/>(with replacement)</FONT>>,
                fillcolor='#F2F3F4', color='#95A5A6', style='filled,dashed']

    seq_choice [shape=diamond, fillcolor='#117864', fontcolor='white', penwidth=2,
                label=<<B>--read-simulator</B>>]

    seq_perfect [label=<<B>perfect</B><BR/><FONT POINT-SIZE='9'>Trim to --read-length (120bp)<BR/>No sequencing errors<BR/>fastq output: dummy Q40</FONT>>,
                 fillcolor='#D1F2EB', color='#117864']

    seq_art [label=<<B>art</B><BR/><FONT POINT-SIZE='9'>ART-modern, Illumina model<BR/>--sequencer-profile (e.g. HiSeq2500_150bp)<BR/>--read-length; optional --paired-end<BR/>(--pe-frag-len-mean/sd)<BR/>renamed via SAM RNAME field</FONT>>,
             fillcolor='#A9DFBF', color='#117864']

    seq_badread [label=<<B>badread</B><BR/><FONT POINT-SIZE='9'>Long-read model: ONT or PacBio CLR<BR/>--sequencer-profile ont / ont-2020 / pacbio<BR/>--coverage-depth (1.0 ~ 1 read/fragment)<BR/>--long-read-length-mean/sd<BR/>--badread-glitches / junk-reads /<BR/>random-reads / chimeras<BR/>renamed via {ref},{start}-{end} in description</FONT>>,
                 fillcolor='#73C6B6', color='#117864']

    seq_merge [label='', shape=point, width=0.15, fillcolor='#117864', color='#117864']

    seq_sample -> seq_choice [penwidth=1.8]
    seq_choice -> seq_perfect [label='perfect', color='#117864']
    seq_choice -> seq_art     [label='art',     color='#117864']
    seq_choice -> seq_badread [label='badread', color='#117864']
    seq_perfect -> seq_merge [color='#117864']
    seq_art     -> seq_merge [color='#117864']
    seq_badread -> seq_merge [color='#117864']
  }

  fs_output:s -> seq_sample:e [penwidth=2, constraint=false]
  td_weight -> seq_sample [style=invis]

  seq_output [label=<<B>reads.fa / reads.fastq</B><BR/><FONT POINT-SIZE='10'>(or reads_R1.fa + reads_R2.fa<BR/>if art --paired-end)</FONT>>,
              fillcolor='#D1F2EB', color='#117864', penwidth=2.5, shape=folder]

  seq_merge -> seq_output [penwidth=2.5, color='#117864']

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 7 — FILTER (optional host depletion)
  // ══════════════════════════════════════════════════════════════════════
  filt_step [label=<<B>baitbench filter</B><BR/><FONT POINT-SIZE='9'>optional; --host genome.fa<BR/>rammap preset sr, removes<BR/>host-aligned reads</FONT>>,
             fillcolor='#B9770E', fontcolor='white', style='filled,dashed', penwidth=1.8]

  seq_output:w -> filt_step:s [style=dashed, penwidth=1.8, constraint=false,
                           label=<<FONT POINT-SIZE='8'>if --host given</FONT>>]
  sim_align -> filt_step [style=invis]

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 8 — MAP + LIST
  // ══════════════════════════════════════════════════════════════════════
  map_step [label=<<B>baitbench map</B><BR/><FONT POINT-SIZE='9'>rammap: reads vs mapping_reference.fa<BR/>preset sr &rarr; mapped.sam<BR/>(paired-end via reads_r2)</FONT>>,
            fillcolor='#154360', fontcolor='white', penwidth=2]

  list_step [label=<<B>baitbench list</B><BR/><FONT POINT-SIZE='9'>Parse SAM &rarr; per-reference<BR/>read counts (detected.list)</FONT>>,
             fillcolor='#0D3349', fontcolor='white', penwidth=2]

  seq_output:w -> map_step:s [penwidth=1.8, constraint=false, label=<<FONT POINT-SIZE='8'>no --host</FONT>>]
  filt_step -> map_step [style=dashed, penwidth=1.8]
  prep_mapping -> map_step [style=dashed, color='#95A5A6']
  map_step -> list_step [penwidth=2]

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 9 — METRICS  (genome-aware 3-way classification)
  // ══════════════════════════════════════════════════════════════════════
  subgraph cluster_metrics {
    label=<<B>baitbench metrics — genome-aware classification</B>>;
    style=rounded; color='#0B5345'; fontsize=13; fontcolor='#0B5345'; margin=16

    met_table [shape=plaintext, label=<
      <TABLE BORDER='1' CELLBORDER='1' CELLSPACING='0' CELLPADDING='5' COLOR='#BDC3C7'>
        <TR><TD BGCOLOR='#D6EAF8'><B>Category</B></TD><TD BGCOLOR='#D6EAF8'><B>Detected</B></TD><TD BGCOLOR='#D6EAF8'><B>Class</B></TD></TR>
        <TR><TD>Sample target</TD><TD>Yes</TD><TD BGCOLOR='#D5F5E3'>TP</TD></TR>
        <TR><TD>Sample target</TD><TD>No</TD><TD BGCOLOR='#FADBD8'>FN</TD></TR>
        <TR><TD>Non-sample target</TD><TD>Yes</TD><TD BGCOLOR='#FADBD8'>FP_target</TD></TR>
        <TR><TD>Non-sample target</TD><TD>No</TD><TD BGCOLOR='#D5F5E3'>TN_target</TD></TR>
        <TR><TD>Distractor</TD><TD>Yes</TD><TD BGCOLOR='#FADBD8'>FP_distractor</TD></TR>
        <TR><TD>Distractor</TD><TD>No</TD><TD BGCOLOR='#D5F5E3'>TN_distractor</TD></TR>
        <TR><TD>Untargeted genome</TD><TD>-</TD><TD BGCOLOR='#F2F3F4'>untargeted</TD></TR>
      </TABLE>
    >]

    met_calc [label=<<B>Sensitivity, Specificity,<BR/>Precision, F1</B><BR/><FONT POINT-SIZE='9'>+ read-level flow:<BR/>sample/nonsample/distractor captured<BR/>correctly vs. incorrectly mapped</FONT>>,
              fillcolor='#0B5345', fontcolor='white', penwidth=2]

    met_out [label=<<B>results.tsv &middot; detected_detail.tsv<BR/>results.json &middot; coverage.tsv</B>>,
             fillcolor='#D5F5E3', color='#27AE60', penwidth=1.8]

    met_table -> met_calc [penwidth=1.8]
    met_calc -> met_out [penwidth=1.8]
  }

  list_step -> met_table [penwidth=2]
  prep_stmap -> met_table [style=dashed, color='#95A5A6']
  prep_idlists -> met_table [style=dashed, color='#95A5A6']

  // ══════════════════════════════════════════════════════════════════════
  // STAGE 10 — REPORT
  // ══════════════════════════════════════════════════════════════════════
  rep_step [label=<<B>baitbench report</B><BR/><FONT POINT-SIZE='9'>Rscript R/report.R + report.Rmd<BR/>Sankey &middot; detection lollipop<BR/>coverage plots &middot; metrics tables<BR/>parameters panel (R / ggplot2)</FONT>>,
            fillcolor='#6C3483', fontcolor='white', penwidth=2]

  rep_output [label=<<B>report.html</B><BR/><FONT POINT-SIZE='10'>Self-contained interactive report</FONT>>,
              fillcolor='#F3E5F5', color='#7D3C98', penwidth=2.5, shape=folder]

  met_out -> rep_step [penwidth=2.5, color='#7D3C98']
  rep_step -> rep_output [penwidth=2.5, color='#7D3C98']
}
")

save_diagram(diagram12, "paper_full_simulation_pipeline.png", width = 6000)
}



message("\nDone. Diagrams saved to: ", outdir)
diagram_names <- c(
  "1. prepare_mode1_standard_nosample.png",
  "2. prepare_mode2_standard_sample.png",
  "3. prepare_mode3_genomes_nosample.png",
  "4. prepare_mode4_genomes_sample.png",
  "5. pipeline_overview.png",
  "6. pipeline_detailed.png",
  "7. paper_workflow_overview.png",
  "8. paper_thermodynamic_scoring.png",
  "9. paper_fragment_sampling.png",
  "10. paper_build_probes.png",
  "11. paper_assess_probes.png",
  "12. paper_full_simulation_pipeline.png"
)
for (i in selected) {
  if (i <= length(diagram_names)) message("  ", diagram_names[i])
}
