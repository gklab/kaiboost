"""Load the proband's HPO terms from a local, uncommitted config file.

The phenotype comes from the hackathon's gated clinical document, which may not
be redistributed, so it must never be hardcoded in a public repository. Put the
terms in config/phenotype.txt (one per line, '#' comments allowed); see
config/phenotype.example.txt for the format.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # mva-analysis
CONFIG = ROOT.parent / "config" / "phenotype.txt"


def load_hpo(path: Path = CONFIG) -> list[str]:
    if not path.exists():
        raise SystemExit(
            f"Phenotype config not found: {path}\n"
            "Copy config/phenotype.example.txt to config/phenotype.txt and fill in\n"
            "the proband's HPO terms from the challenge's clinical document.\n"
            "This file is gitignored on purpose - do not commit patient phenotype."
        )
    terms = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            if not line.startswith("HP:"):
                raise SystemExit(f"{path}: not an HPO id: {line!r}")
            terms.append(line)
    if not terms:
        raise SystemExit(f"{path}: no HPO terms found")
    return terms
