#!/bin/bash
# Negative control: GIAB HG002 (karyotypically normal male, PCR-free TruSeq,
# 2x148) through the IDENTICAL pipeline as the proband. If HG002 reproduces the
# proband's chr16/17/18/19/21 elevation, that pattern is a pipeline artifact;
# if HG002 comes out flat, the proband's pattern needs a biological explanation.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$ROOT/control_pass"
MMI="$ROOT/cache/genome.sr.mmi"
BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1"
PART=$((250*1024*1024))   # 250 MB per mate per lane, matching the proband pass
mkdir -p "$WORK"; cd "$WORK"

fetch_part() {  # url, out, bytes
  local url=$1 out=$2 part=$3 have
  while :; do
    have=$(stat -f%z "$out" 2>/dev/null || echo 0)
    [ "$have" -ge "$part" ] && break
    curl -sfL --connect-timeout 20 --speed-time 60 --speed-limit 5000 \
      -r "$have-$((part-1))" "$url" >> "$out" || sleep 3
  done
}

for lane in L001 L002 L003 L004; do
  [ -f "ctrl_$lane.bins.tsv" ] && { echo "[ctrl] $lane done"; continue; }
  echo "[ctrl] $lane downloading $(date +%H:%M)"
  fetch_part "$BASE/2A1_CGATGT_${lane}_R1_001.fastq.gz" "$lane.R1.gz" "$PART"
  fetch_part "$BASE/2A1_CGATGT_${lane}_R2_001.fastq.gz" "$lane.R2.gz" "$PART"
  n1=$( (gunzip -c "$lane.R1.gz" 2>/dev/null || true) | wc -l )
  n2=$( (gunzip -c "$lane.R2.gz" 2>/dev/null || true) | wc -l )
  n=$(( (n1 < n2 ? n1 : n2) / 4 * 4 ))
  echo "[ctrl] $lane aligning $((n/4)) pairs $(date +%H:%M)"
  minimap2 -ax sr -t 11 --secondary=no "$MMI" \
    <( (gunzip -c "$lane.R1.gz" 2>/dev/null || true) | head -n "$n" ) \
    <( (gunzip -c "$lane.R2.gz" 2>/dev/null || true) | head -n "$n" ) \
    2> "$lane.mm2.log" | python3 "$ROOT/tools/stream_filter.py" "ctrl_$lane"
  [ -s "ctrl_$lane.bins.tsv" ] && rm -f "$lane.R1.gz" "$lane.R2.gz"
done
echo "CONTROL_DONE $(date)"
