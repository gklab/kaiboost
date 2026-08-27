"""Shared helpers: paths, disk cache, Ensembl / MyVariant API clients, bcftools wrapper."""
import json
import hashlib
import subprocess
import time
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT / "cache"
RESULTS_DIR = ROOT / "results"
VCF = ROOT.parent / "mva-data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
PROBAND_ID = "WGS_EX2312012"

ENSEMBL = "https://rest.ensembl.org"
MYVARIANT = "https://myvariant.info/v1"


def _cache_path(namespace: str, key: str) -> Path:
    h = hashlib.sha256(key.encode()).hexdigest()[:24]
    d = CACHE_DIR / namespace
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{h}.json"


def cached_json(namespace: str, key: str, fetch, ttl_days: int = 365):
    """Fetch-through disk cache. `fetch` is a zero-arg callable returning a JSON-able object."""
    p = _cache_path(namespace, key)
    if p.exists() and (time.time() - p.stat().st_mtime) < ttl_days * 86400:
        return json.loads(p.read_text())
    val = fetch()
    p.write_text(json.dumps(val))
    return val


def http_json(url: str, data: dict | None = None, timeout: int = 120, retries: int = 3):
    """GET (data=None) or POST JSON with retries/backoff."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                data=json.dumps(data).encode() if data is not None else None,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            code = getattr(e, "code", None)
            if code == 404:
                return None
            if code == 400:  # bad payload — retrying is pointless; surface the server's message
                body = ""
                try:
                    body = e.read().decode()[:500]
                except Exception:
                    pass
                raise ValueError(f"HTTP 400 from {url.split('?')[0]}: {body}") from e
            if attempt == retries - 1:
                raise
            wait = 30 if code == 429 else 2 ** attempt * 2
            time.sleep(wait)


def gene_coords(symbol: str) -> dict | None:
    """GRCh38 coordinates for a gene symbol via Ensembl (cached). Returns dict or None."""
    def fetch():
        r = http_json(f"{ENSEMBL}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json")
        if r is None:
            return {}
        return {
            "symbol": r.get("display_name", symbol),
            "chrom": r["seq_region_name"],
            "start": r["start"],
            "end": r["end"],
            "strand": r["strand"],
            "description": r.get("description", ""),
        }
    r = cached_json("gene_coords", symbol.upper(), fetch)
    return r or None


def vep_region(variants: list[str], params: str = "af=1&af_gnomadg=1&hgvs=1&numbers=1&canonical=1"):
    """Annotate VCF-style variant strings ('15 40209701 . T G . . .') via VEP REST (cached).
    Handles batching (max 200/request)."""
    out = []
    for i in range(0, len(variants), 200):
        chunk = variants[i:i + 200]
        key = params + "|" + "|".join(chunk)
        out.extend(cached_json("vep", key,
                               lambda c=chunk: http_json(f"{ENSEMBL}/vep/homo_sapiens/region?{params}",
                                                         {"variants": c}) or []))
    return out


def myvariant(hgvs_ids: list[str]):
    """Batch query MyVariant.info (hg38) for gnomAD/ClinVar/CADD/REVEL etc. (cached).
    hgvs_ids like 'chr15:g.40209701T>G'. Max 1000/request."""
    fields = ("gnomad_genome.af,gnomad_exome.af,clinvar.rcv.clinical_significance,"
              "clinvar.clinical_significance,cadd.phred,dbnsfp.revel.score,"
              "dbnsfp.alphamissense.score,dbnsfp.spliceai.ds_max,dbsnp.rsid,"
              "dbnsfp.genename,snpeff.ann")
    out = []
    for i in range(0, len(hgvs_ids), 1000):
        chunk = hgvs_ids[i:i + 1000]
        key = fields + "|" + "|".join(chunk)

        def fetch(c=chunk):
            body = ("ids=" + ",".join(c) + f"&fields={fields}&assembly=hg38")
            for attempt in range(4):
                try:
                    req = urllib.request.Request(f"{MYVARIANT}/variant", data=body.encode(),
                                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
                    with urllib.request.urlopen(req, timeout=180) as r:
                        return json.load(r)
                except urllib.error.HTTPError as e:
                    if e.code == 429 and attempt < 3:
                        time.sleep(30 * (attempt + 1))
                        continue
                    raise
        out.extend(cached_json("myvariant", key, fetch))
    return out


def bcftools(args: list[str]) -> str:
    """Run bcftools, return stdout."""
    r = subprocess.run(["bcftools"] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"bcftools {' '.join(args)} failed:\n{r.stderr}")
    return r.stdout


def vcf_variants_in_region(region: str, pass_only: bool = False) -> list[dict]:
    """Extract variants in a region (e.g. '15:40160984-40221137'). Multiallelics are split."""
    args = ["query", "-r", region,
            "-f", r"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t[%GT]\t[%AD]\t[%DP]\t[%GQ]\n"]
    if pass_only:
        args += ["-i", 'FILTER="PASS"']
    args.append(str(VCF))
    rows = []
    for line in bcftools(args).strip().split("\n"):
        if not line:
            continue
        c, pos, vid, ref, alts, filt, gt, ad, dp, gq = line.split("\t")
        for alt in alts.split(","):
            if alt in ("<NON_REF>", "*"):
                continue
            rows.append(dict(chrom=c, pos=int(pos), id=vid, ref=ref, alt=alt,
                             filter=filt, gt=gt, ad=ad, dp=dp, gq=gq))
    return rows
