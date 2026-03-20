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

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Parse diagram selection (e.g. "5", "2-4", "1,3,5-6")
parse_selection <- function(spec, total = 6) {
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

save_diagram <- function(graph, filename, width = 3200) {
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


message("\nDone. Diagrams saved to: ", outdir)
diagram_names <- c(
  "1. prepare_mode1_standard_nosample.png",
  "2. prepare_mode2_standard_sample.png",
  "3. prepare_mode3_genomes_nosample.png",
  "4. prepare_mode4_genomes_sample.png",
  "5. pipeline_overview.png",
  "6. pipeline_detailed.png"
)
for (i in selected) {
  if (i <= length(diagram_names)) message("  ", diagram_names[i])
}
