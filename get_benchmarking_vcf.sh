#!/usr/bin/env bash
# =============================================================================
# wgs_to_wes_filter.sh
#
# Simulates a WES-like callset from the NA12878 (HG001) WGS benchmark VCF.
#
# Strategy to approximate real WES variant counts (~30-40K):
#   1. Restrict to CDS regions only (not full exons) — excludes UTRs
#   2. Restrict to protein-coding genes only — excludes lncRNA, pseudogenes etc.
#   3. Merge overlapping intervals before filtering
#
# This brings the target BED from ~90 Mb (all exons) to ~35 Mb (coding only),
# which closely matches commercial WES capture kits (Agilent, Illumina, Twist).
#
# Dependencies: wget, bcftools, bedtools, bgzip, tabix, awk, sort, gzip
# =============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Colour helpers for readable log output
# ------------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[1;33m'; NC='\033[0m'
log()  { echo -e "${CYAN}[$(date '+%H:%M:%S')] $*${NC}"; }
ok()   { echo -e "${GREEN}[$(date '+%H:%M:%S')] $*${NC}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARNING: $*${NC}"; }
err()  { echo -e "${RED}[ERROR] $*${NC}" >&2; exit 1; }

# ------------------------------------------------------------------------------
# Source URLs
# ------------------------------------------------------------------------------

# GIAB NA12878 (HG001) benchmark VCF — GRCh38, chromosomes 1-22
VCF_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TBI_URL="${VCF_URL}.tbi"

# GENCODE v46 comprehensive annotation (GRCh38)
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"

# ------------------------------------------------------------------------------
# File names
# ------------------------------------------------------------------------------
VCF_FILE="HG001_GRCh38.vcf.gz"          # Downloaded WGS VCF
TBI_FILE="${VCF_FILE}.tbi"               # Its tabix index
GTF_FILE="gencode.v46.gtf.gz"           # Downloaded GTF annotation
RAW_BED="cds_protein_coding_raw.bed"    # CDS intervals before merging
MERGED_BED="cds_protein_coding_merged.bed"      # After bedtools merge
MERGED_BED_GZ="${MERGED_BED}.gz"        # bgzipped + tabix-indexed version
OUT_VCF="HG001.WES_like.vcf.gz"         # Final filtered output

# ------------------------------------------------------------------------------
# Dependency check — fail early if anything is missing
# ------------------------------------------------------------------------------
log "Checking dependencies..."
for cmd in wget bcftools bedtools bgzip tabix awk sort gzip; do
    command -v "$cmd" >/dev/null 2>&1 \
        || err "'$cmd' not found. Please install it and re-run."
done
ok "All dependencies found."

# ==============================================================================
# STEP 1 — Download WGS VCF + tabix index
# ==============================================================================
log "[1/4] Downloading NA12878 benchmark VCF (GRCh38)..."
wget -q --show-progress -O "$VCF_FILE" "$VCF_URL"
wget -q --show-progress -O "$TBI_FILE" "$TBI_URL"
ok "VCF downloaded -> ${VCF_FILE}"

# Count total variants before filtering (used for retention rate at the end)
log "Counting variants in source WGS VCF..."
TOTAL_VARIANTS=$(bcftools view -H "$VCF_FILE" | wc -l)
log "Total WGS variants: ${TOTAL_VARIANTS}"

# ==============================================================================
# STEP 2 — Download GENCODE GTF annotation
# ==============================================================================
log "[2/4] Downloading GENCODE v46 GTF annotation..."
wget -q --show-progress -O "$GTF_FILE" "$GTF_URL"
ok "GTF downloaded -> ${GTF_FILE}"

# ==============================================================================
# STEP 3 — Build a CDS-only, protein-coding BED file
# ==============================================================================
log "[3/4] Building CDS BED for protein-coding genes..."

# Why CDS instead of "exon"?
#   GTF "exon" features include 5' and 3' UTRs, which real WES kits mostly
#   don't capture. "CDS" (coding sequence) covers only the translated region,
#   matching the ~35 Mb targeted by commercial kits like Agilent SureSelect
#   or Twist Human Core Exome.
#
# Why protein_coding only?
#   GENCODE annotates ~60K transcripts including lncRNAs, pseudogenes, and
#   other non-coding RNAs. Restricting to protein_coding genes drops the
#   target space from ~90 Mb to ~35 Mb, consistent with real WES.
#
# Coordinate conversion:
#   GTF is 1-based inclusive [start, end]
#   BED is 0-based half-open  [start-1, end)

log "  Extracting CDS intervals from protein-coding genes..."
gzip -dc "$GTF_FILE" \
  | awk '
      # Select only CDS features on canonical chromosomes
      $3 == "CDS" && $1 ~ /^chr([0-9]{1,2}|X|Y|M)$/ {

          # Only keep protein-coding genes
          # The gene_type attribute is in the last (9th) field of the GTF
          if ($0 ~ /gene_type "protein_coding"/) {

              # Convert GTF 1-based start to BED 0-based start ($4 - 1)
              # End coordinate is the same in both formats ($5)
              print $1 "\t" ($4 - 1) "\t" $5
          }
      }
    ' \
  | sort -k1,1V -k2,2n \
  > "$RAW_BED"

RAW_INTERVALS=$(wc -l < "$RAW_BED")
log "  Raw CDS intervals (pre-merge): ${RAW_INTERVALS}"

# Merge overlapping/adjacent intervals.
# bcftools --regions-file requires non-overlapping intervals; overlaps cause
# duplicate variant records. bedtools merge collapses them correctly.
log "  Merging overlapping CDS intervals with bedtools..."
bedtools merge -i "$RAW_BED" > "$MERGED_BED"

MERGED_INTERVALS=$(wc -l < "$MERGED_BED")
log "  Merged CDS intervals: ${MERGED_INTERVALS}"

# Calculate approximate target size in Mb for sanity check
TARGET_MB=$(awk '{sum += $3 - $2} END {printf "%.1f", sum/1e6}' "$MERGED_BED")
log "  Total target size: ${TARGET_MB} Mb  (real WES kits: ~33-65 Mb)"

# bgzip-compress and tabix-index the BED.
# bcftools --regions-file requires a tabix-indexed, bgzipped BED for
# random access; plain text BED files are not accepted with this flag.
log "  Compressing and indexing BED..."
bgzip -f "$MERGED_BED"
tabix -p bed "$MERGED_BED_GZ"
ok "CDS BED ready -> ${MERGED_BED_GZ}  (${TARGET_MB} Mb, ${MERGED_INTERVALS} intervals)"

# ==============================================================================
# STEP 4 — Filter WGS VCF to CDS regions
# ==============================================================================
log "[4/4] Filtering WGS VCF to protein-coding CDS regions..."

# --regions-file: uses the tabix index for efficient random-access lookup.
#   This is faster than --targets-file (which streams the whole VCF).
# -Oz: output as bgzipped VCF
bcftools view \
    --regions-file "$MERGED_BED_GZ" \
    --output-type z \
    --output "$OUT_VCF" \
    "$VCF_FILE"

# Index the output for downstream tools (GATK, IGV, etc.)
bcftools index --tbi "$OUT_VCF"
ok "Filtered VCF written -> ${OUT_VCF}"

# ==============================================================================
# SUMMARY
# ==============================================================================
WES_VARIANTS=$(bcftools view -H "$OUT_VCF" | wc -l)
PCTG=$(awk "BEGIN{printf \"%.1f\", 100*${WES_VARIANTS}/${TOTAL_VARIANTS}}")

echo ""
echo "============================================================"
echo "  WGS to WES-like filtering — COMPLETE"
echo "============================================================"
printf "  %-30s %s\n"   "Source WGS VCF:"          "$VCF_FILE"
printf "  %-30s %s\n"   "CDS BED (bgzipped):"      "$MERGED_BED_GZ"
printf "  %-30s %s\n"   "Output WES-like VCF:"     "$OUT_VCF"
echo "  ------------------------------------------------------------"
printf "  %-30s %s Mb\n" "CDS target size:"        "$TARGET_MB"
printf "  %-30s %d\n"   "Merged CDS intervals:"    "$MERGED_INTERVALS"
echo "  ------------------------------------------------------------"
printf "  %-30s %d\n"   "Variants (WGS total):"    "$TOTAL_VARIANTS"
printf "  %-30s %d\n"   "Variants (WES-like):"     "$WES_VARIANTS"
printf "  %-30s %s%%\n"  "Retention rate:"          "$PCTG"
echo "============================================================"
