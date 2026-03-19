#!/usr/bin/env Rscript
# BaitBench Prepare Step Diagrams
# Generates 4 diagrams showing inputs, outputs, and fragment composition
# for each combination of --sample and --genomes flags.
#
# Usage: Rscript R/prepare_diagrams.R [output_dir]
# Output: 4 PNG files in the output directory (wide format for slides)

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save_diagram <- function(graph, filename) {
  svg <- export_svg(graph)
  rsvg_png(charToRaw(svg), file = file.path(outdir, filename), width = 3200)
  message("Saved: ", file.path(outdir, filename))
}


# ============================================================================
# Diagram 1: Standard Mode, No --sample
# ============================================================================
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


# ============================================================================
# Diagram 2: Standard Mode, With --sample
# ============================================================================
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


# ============================================================================
# Diagram 3: Genome Mode, No --sample
# ============================================================================
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


# ============================================================================
# Diagram 4: Genome Mode, With --sample
# ============================================================================
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


message("\nAll 4 diagrams saved to: ", outdir)
message("  1. prepare_mode1_standard_nosample.png")
message("  2. prepare_mode2_standard_sample.png")
message("  3. prepare_mode3_genomes_nosample.png")
message("  4. prepare_mode4_genomes_sample.png")
