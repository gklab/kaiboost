#!/bin/bash
# Negative control, part 2 — adds units until the HG002 control matches the
# proband's depth (~13M read pairs). Units are independent gzip streams
# (separate flowcells / file parts), so each can be prefix-truncated safely.
#
# Fixes over control_pass.sh: HEAD-check every URL before fetching (flowcell
# AH8VC6ADXX is a 2-lane rapid-run card and has no L003/L004 — the old loop
# spun forever on the 404), and cap retries so a dead unit is skipped, not hung.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$ROOT/control_pass"
MMI="$ROOT/cache/genome.sr.mmi"
R="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq"
S="Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT"
mkdir -p "$WORK"; cd "$WORK"

# label|flowcell|lane|part|max_MB   (ctrl_L001/ctrl_L002 already done: 0018 _001)
UNITS="
ctrl_u3|140528_D00360_0018_AH8VC6ADXX|L001|002|999
ctrl_u4|140528_D00360_0018_AH8VC6ADXX|L002|002|999
ctrl_u5|140528_D00360_0019_BH8VDAADXX|L001|001|250
ctrl_u6|140528_D00360_0019_BH8VDAADXX|L002|001|250
"

remote_size() { curl -sIL --connect-timeout 15 "$1" | tr -d '\r' \
  | awk 'tolower($1)=="content-length:"{v=$2} END{print v+0}'; }

fetch_part() {  # url, out, bytes -> 0 ok / 1 gave up
  local url=$1 out=$2 part=$3 have tries=0
  while :; do
    have=$(stat -f%z "$out" 2>/dev/null || echo 0)
    [ "$have" -ge "$part" ] && return 0
    tries=$((tries+1))
    [ "$tries" -gt 40 ] && { echo "  [give up after 40 tries] $url"; return 1; }
    curl -sfL --connect-timeout 20 --speed-time 60 --speed-limit 5000 \
      -r "$have-$((part-1))" "$url" >> "$out" || sleep 3
  done
}

for u in $UNITS; do
  IFS='|' read -r lbl fc lane part maxmb <<< "$u"
  [ -f "$lbl.bins.tsv" ] && { echo "[ctrl] $lbl already done"; continue; }
  u1="$R/$fc/${S}_${lane}_R1_${part}.fastq.gz"
  u2="$R/$fc/${S}_${lane}_R2_${part}.fastq.gz"
  s1=$(remote_size "$u1"); s2=$(remote_size "$u2")
  if [ "${s1:-0}" -lt 1000000 ] || [ "${s2:-0}" -lt 1000000 ]; then
    echo "[ctrl] $lbl SKIP — missing on server (R1=$s1 R2=$s2)"; continue
  fi
  cap=$((maxmb*1024*1024))
  p1=$(( s1 < cap ? s1 : cap )); p2=$(( s2 < cap ? s2 : cap ))
  echo "[ctrl] $lbl downloading $((p1/1048576))+$((p2/1048576))MB $(date +%H:%M)"
  fetch_part "$u1" "$lbl.R1.gz" "$p1" || { echo "[ctrl] $lbl download failed"; continue; }
  fetch_part "$u2" "$lbl.R2.gz" "$p2" || { echo "[ctrl] $lbl download failed"; continue; }
  n1=$( (gunzip -c "$lbl.R1.gz" 2>/dev/null || true) | wc -l )
  n2=$( (gunzip -c "$lbl.R2.gz" 2>/dev/null || true) | wc -l )
  n=$(( (n1 < n2 ? n1 : n2) / 4 * 4 ))
  echo "[ctrl] $lbl aligning $((n/4)) pairs $(date +%H:%M)"
  minimap2 -ax sr -t 11 --secondary=no "$MMI" \
    <( (gunzip -c "$lbl.R1.gz" 2>/dev/null || true) | head -n "$n" ) \
    <( (gunzip -c "$lbl.R2.gz" 2>/dev/null || true) | head -n "$n" ) \
    2> "$lbl.mm2.log" | python3 "$ROOT/tools/stream_filter.py" "$lbl"
  [ -s "$lbl.bins.tsv" ] && rm -f "$lbl.R1.gz" "$lbl.R2.gz"
done
rm -f L003.R1.gz L003.R2.gz L004.R1.gz L004.R2.gz
echo "CONTROL2_DONE $(date)"
