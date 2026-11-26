#!/usr/bin/env bash
# Fetch the D6331 mock-community reads from ENA.
#
# Why: every SimVal comparison so far has been against pre-computed *count tables*, so
# the real side was mapped by whatever tool produced them and the simulated side by
# BaitBench's embedded aligner. Any mapper difference sits inside the residual,
# unattributed. With reads in hand both sides can be mapped identically.
#
# ENA rather than SRA-tools: it serves the FASTQs directly over HTTPS, so there is
# nothing to install and nothing to convert.
#
# Download paths come from ENA's filereport API, not from a rule applied to the
# accession. The layout is `<first 6 chars>/0<last 2 digits>/<accession>/`, which is
# easy to get wrong — reconstructing it as the last *three* digits 404s on every run.
#
# NOTE on the filenames: ENA serves these as `<acc>_subreads.fastq.gz`, but the content
# is **CCS/HiFi**, not subreads — the read headers end `.../ccs/N`. SRA applies the
# `_subreads` suffix generically to PacBio submissions. Checked rather than assumed,
# because the distinction decides whether the read length tracks the insert directly
# (CCS) or is biased toward short inserts by pass count (subreads).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config.sh"

# SIMVAL_READS defaults (config.sh) to the directory these reads already live in,
# BB_Paper/Comp/sra/fastq — they were fetched once for map_real_reads.sh. This script
# exists to re-fetch them if that copy is ever lost, and to record the accessions; it
# skips anything already present, so running it on a complete directory is a no-op.
# Files there are named <sample>_<accession>.fastq.gz; a bare <sample>.fastq.gz written
# by this script is treated as equivalent.
DEST="${SIMVAL_READS:-$SIMVAL_DIR/reads}"
PROJECT="${SIMVAL_BIOPROJECT:-PRJNA1081843}"
mkdir -p "$DEST"

acc_file="$HERE/accessions.tsv"
[[ -s "$acc_file" ]] || { echo "missing $acc_file" >&2; exit 1; }

report="$DEST/.ena_filereport.tsv"
if [[ ! -s "$report" ]]; then
    echo "querying ENA for $PROJECT ..."
    curl -fsSL --retry 3 \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,read_count&format=tsv" \
        -o "$report"
fi

fetched=0; skipped=0; failed=0
while IFS=$'\t' read -r sample acc; do
    [[ -z "${sample:-}" || "$sample" == \#* || "$sample" == "sample" ]] && continue
    out="$DEST/${sample}.fastq.gz"
    # Either naming counts as present.
    if [[ -s "$out" ]] || compgen -G "$DEST/${sample}_*.fastq.gz" > /dev/null; then
        skipped=$((skipped+1)); continue
    fi

    path=$(awk -F'\t' -v a="$acc" '$1==a {print $2}' "$report" | cut -d';' -f1)
    if [[ -z "$path" ]]; then
        echo "  no ENA path for $acc ($sample)" >&2; failed=$((failed+1)); continue
    fi
    echo "[$(date +%H:%M:%S)] $sample <- $acc"
    # ENA reports the path with its FTP host prefix; swap that for HTTPS on the same host.
    if ! curl -fsSL --retry 3 --retry-delay 5 -o "$out.part" \
            "https://ftp.sra.ebi.ac.uk/${path#ftp.sra.ebi.ac.uk/}"; then
        echo "  FAILED: $path" >&2; rm -f "$out.part"; failed=$((failed+1)); continue
    fi
    mv "$out.part" "$out"
    fetched=$((fetched+1))
done < "$acc_file"

echo
echo "fetched $fetched, already present $skipped, failed $failed -> $DEST"
du -sh "$DEST" 2>/dev/null
