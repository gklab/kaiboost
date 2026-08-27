#!/bin/bash
# 10% uniform sample pass over all 4 lanes: resumable range downloads (prefetch
# next lane while aligning current), minimap2 full-genome, streaming filter.
# Products per lane: fastq_pass/<lane>.bins.tsv + <lane>.locus.sam
# Run under: caffeinate -is bash tools/sample_pass.sh
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CACHE="$ROOT/cache"
WORK="$ROOT/fastq_pass"
MMI="$CACHE/genome.sr.mmi"
FRAC_NUM=40   # take first 1/40 (2.5%) of each file — chromosome-level mosaic
              # sensitivity is GC-systematics-bound, not read-count-bound
REPO="SageBio/mva-hackathon-2026-data"
LANES="L001 L002 L003 L004"
mkdir -p "$WORK"
cd "$WORK"

TOK=$(python3 -c "from huggingface_hub import get_token; print(get_token())")

url_for() {  # lane, mate
  python3 -c "from huggingface_hub import hf_hub_url; print(hf_hub_url('$REPO', 'WGS_EX2312012_HGWCNDSX7_S16_${1}_R${2}_001.fastq.gz', repo_type='dataset'))"
}

size_for() {  # url -> content-length
  curl -sIL -H "Authorization: Bearer $TOK" "$1" | tr -d '\r' \
    | awk 'tolower($1)=="content-length:"{s=$2} END{print s}'
}

fetch_part() {  # url, out, part_bytes : resumable range download loop
  local url=$1 out=$2 part=$3 have
  while :; do
    have=$(stat -f%z "$out" 2>/dev/null || echo 0)
    [ "$have" -ge "$part" ] && break
    curl -sfL --connect-timeout 20 --speed-time 60 --speed-limit 5000 \
      -r "$have-$((part-1))" -H "Authorization: Bearer $TOK" "$url" >> "$out"
    sleep 3
  done
}

download_lane() {  # lane : fetch R1+R2 parts concurrently, mark .done
  local lane=$1 u1 u2 s1 s2
  [ -f "$lane.done" ] && return
  u1=$(url_for "$lane" 1); u2=$(url_for "$lane" 2)
  s1=$(size_for "$u1");    s2=$(size_for "$u2")
  echo "[dl] $lane R1 $((s1/FRAC_NUM/1048576))MB R2 $((s2/FRAC_NUM/1048576))MB"
  fetch_part "$u1" "$lane.R1.part.gz" "$((s1/FRAC_NUM))" &
  fetch_part "$u2" "$lane.R2.part.gz" "$((s2/FRAC_NUM))" &
  wait
  touch "$lane.done"
  echo "[dl] $lane complete"
}

# background prefetcher: lanes strictly in order (bandwidth is one shared pipe)
( for lane in $LANES; do download_lane "$lane"; done ) &
PREFETCH=$!

for lane in $LANES; do
  until [ -f "$lane.done" ]; do sleep 20; done
  if [ -f "$lane.bins.tsv" ]; then echo "[align] $lane already done"; continue; fi
  echo "[align] $lane starting $(date +%H:%M)"
  # truncated gz tails are expected; pair-sync via equal head line counts
  n1=$( (gunzip -c "$lane.R1.part.gz" 2>/dev/null || true) | wc -l )
  n2=$( (gunzip -c "$lane.R2.part.gz" 2>/dev/null || true) | wc -l )
  nreads=$(( (n1 < n2 ? n1 : n2) / 4 * 4 ))
  echo "[align] $lane using $((nreads/4)) read pairs"
  minimap2 -ax sr -t 11 --secondary=no "$MMI" \
    <( (gunzip -c "$lane.R1.part.gz" 2>/dev/null || true) | head -n "$nreads" ) \
    <( (gunzip -c "$lane.R2.part.gz" 2>/dev/null || true) | head -n "$nreads" ) \
    2> "$lane.mm2.log" \
    | python3 "$ROOT/tools/stream_filter.py" "$lane"
  if [ -s "$lane.bins.tsv" ]; then
    rm -f "$lane.R1.part.gz" "$lane.R2.part.gz"
    echo "[align] $lane finished $(date +%H:%M)"
  else
    echo "[align] $lane FAILED (bins missing) — parts kept for retry"
  fi
done
wait "$PREFETCH" 2>/dev/null
echo "ALL_LANES_DONE $(date)"
