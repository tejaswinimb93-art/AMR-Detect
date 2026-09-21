"""
AMRIVA AMRFinderPlus result organizer
-------------------------------------
Drop-in parsing/organization layer for an existing AMRFinderPlus runner.

Use this module AFTER AMRFinderPlus returns its TSV. It:
1. Preserves the technical AMRFinderPlus columns.
2. Normalizes column names.
3. Extracts organism/sample information from optional metadata or FASTA headers.
4. Adds an "Associated antibiotic/group" column from AMRFinderPlus Class/Subclass.
5. Supports cross-sample filtering by antibiotic/group.
6. Keeps genotype wording scientifically cautious: a genomic determinant is not by itself a
   phenotypic AST result.

Expected AMRFinderPlus TSV columns include fields such as:
Element symbol, Element name, Class, Subclass, Method,
% Identity to reference, % Coverage of reference, Contig id, Start, Stop,
Closest reference accession, Closest reference name, etc.
"""

from __future__ import annotations
import re
from pathlib import Path
import pandas as pd


def _clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Map common AMRFinderPlus spelling variants to stable internal names."""
    df = df.copy()
    rename = {}
    aliases = {
        "sample": ["sample", "sample id", "sample_id"],
        "organism": ["organism", "scientific name", "species"],
        "determinant": ["element symbol", "amr determinant", "amr_determinant", "gene"],
        "element_name": ["element name", "mechanism / product", "mechanism/product"],
        "class": ["class"],
        "subclass": ["subclass"],
        "method": ["method"],
        "identity": ["% identity to reference", "identity (%)", "identity"],
        "coverage": ["% coverage of reference", "coverage (%)", "coverage"],
        "contig": ["contig id", "sequence / contig", "sequence/contig", "contig"],
        "start": ["start"],
        "stop": ["stop"],
        "reference": [
            "closest reference accession",
            "reference",
            "reference accession",
        ],
        "raw_result": ["raw result", "raw_result"],
    }
    normalized = {re.sub(r"\s+", " ", str(c).strip().lower()): c for c in df.columns}
    for target, names in aliases.items():
        for name in names:
            key = re.sub(r"\s+", " ", name.lower())
            if key in normalized:
                rename[normalized[key]] = target
                break
    df = df.rename(columns=rename)
    return df


# This mapping is intentionally conservative. "Associated antimicrobial group" is
# preferred to claiming that a gene proves phenotypic resistance to one exact drug.
SUBCLASS_GROUPS = {
    "CARBAPENEM": "Carbapenems",
    "CEPHALOSPORIN": "Cephalosporins",
    "CEPHALOTHIN": "Cephalothin / early cephalosporin group",
    "METHICILLIN": "Methicillin / anti-staphylococcal beta-lactam group",
    "PENICILLIN": "Penicillins",
    "ERYTHROMYCIN": "Erythromycin / macrolide group",
    "TELITHROMYCIN": "Telithromycin / macrolide group",
    "TYLOSIN": "Tylosin / macrolide group",
    "VANCOMYCIN": "Vancomycin / glycopeptide group",
    "CHLORAMPHENICOL": "Chloramphenicol",
    "FLORFENICOL": "Florfenicol",
    "LINEZOLID": "Linezolid / oxazolidinone group",
    "FOSFOMYCIN": "Fosfomycin",
    "COLISTIN": "Colistin",
    "RIFAMPIN": "Rifampin / rifamycin group",
    "RIFAMYCIN": "Rifamycin group",
    "SULFONAMIDE": "Sulfonamides",
    "TETRACYCLINE": "Tetracyclines",
    "TIGECYCLINE": "Tigecycline / tetracycline group",
    "TRIMETHOPRIM": "Trimethoprim",
    "LINCOSAMIDE": "Lincosamides",
    "FLUOROQUINOLONE": "Fluoroquinolones",
    "QUINOLONE": "Quinolones",
}


CLASS_GROUPS = {
    "BETA-LACTAM": "Beta-lactams",
    "MACROLIDE": "Macrolides",
    "GLYCOPEPTIDE": "Glycopeptides",
    "TETRACYCLINE": "Tetracyclines",
    "QUINOLONE": "Quinolones",
    "FLUOROQUINOLONE": "Fluoroquinolones",
    "AMINOGLYCOSIDE": "Aminoglycosides",
    "PHENICOL": "Phenicol antibiotics",
    "SULFONAMIDE": "Sulfonamides",
    "TRIMETHOPRIM": "Trimethoprim",
    "FOSFOMYCIN": "Fosfomycin",
    "COLISTIN": "Colistin",
}


def associated_antibiotic(row) -> str:
    """Return the most specific antimicrobial group supported by Class/Subclass."""
    subclass = _clean(row.get("subclass", "")).upper()
    cls = _clean(row.get("class", "")).upper()

    determinant = _clean(row.get("determinant", "")).lower()

    # Determinant-specific display is used only where the association is well established.
    # It is still labelled as an association, not as a phenotypic AST result.
    if determinant in {"blatem-1", "blaTEM-1".lower(), "blatem"}:
        return "Penicillins; early/narrow-spectrum cephalosporins (associated)"

    if subclass in SUBCLASS_GROUPS:
        return SUBCLASS_GROUPS[subclass]

    if subclass and subclass not in {"AMR", "STRESS", "VIRULENCE"}:
        return subclass.title()

    if cls in CLASS_GROUPS:
        return CLASS_GROUPS[cls]

    return "Not specified by AMRFinderPlus"


def parse_fasta_metadata(fasta_file) -> dict:
    """
    Reads sample/organism from FASTA headers.

    Supported header:
      >S01|organism=Escherichia_coli|...
    or
      >S01 organism=Escherichia_coli ...

    IMPORTANT: a plain nucleotide FASTA does not inherently guarantee an organism name.
    The organism must come from metadata/header or a separate taxonomic identification step.
    """
    path = Path(fasta_file)
    mapping = {}
    current = None
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                sample_id = header.split()[0].split("|")[0]
                m = re.search(r"(?:^|\|)\s*organism=([^|]+)", header, re.I)
                organism = m.group(1).replace("_", " ") if m else ""
                mapping[sample_id] = organism
                current = sample_id
    return mapping


def attach_sample_metadata(df: pd.DataFrame, fasta_file=None, metadata_file=None):
    df = df.copy()

    if "sample" not in df.columns:
        # If AMRFinderPlus was run one file at a time, the caller should pass sample_id.
        df["sample"] = ""

    organism_map = {}

    if metadata_file:
        meta = pd.read_csv(metadata_file)
        cols = {str(c).strip().lower(): c for c in meta.columns}
        sid = cols.get("sample_id") or cols.get("sample")
        org = cols.get("organism") or cols.get("scientific name") or cols.get("species")
        if sid and org:
            organism_map.update(
                dict(zip(meta[sid].astype(str).str.strip(), meta[org].astype(str).str.strip()))
            )

    if fasta_file:
        organism_map.update(parse_fasta_metadata(fasta_file))

    if "organism" not in df.columns:
        df["organism"] = ""

    if organism_map:
        df["organism"] = df.apply(
            lambda r: organism_map.get(_clean(r.get("sample")), _clean(r.get("organism", ""))),
            axis=1,
        )

    df["organism"] = df["organism"].replace("", "Not provided")
    return df


def organize_amrfinder(tsv_file, fasta_file=None, metadata_file=None, sample_id=None):
    """Read and organize one AMRFinderPlus TSV file."""
    df = pd.read_csv(tsv_file, sep="\t", dtype=str, keep_default_na=False)
    df = normalize_columns(df)

    if sample_id:
        df["sample"] = sample_id

    df = attach_sample_metadata(df, fasta_file=fasta_file, metadata_file=metadata_file)

    if "determinant" not in df.columns:
        df["determinant"] = ""
    if "class" not in df.columns:
        df["class"] = ""
    if "subclass" not in df.columns:
        df["subclass"] = ""

    df["Associated antibiotic/group"] = df.apply(associated_antibiotic, axis=1)

    # Keep the useful columns first, then retain every original AMRFinderPlus field.
    preferred = [
        "sample", "organism", "determinant", "Associated antibiotic/group",
        "class", "subclass", "element_name", "method", "identity", "coverage",
        "contig", "start", "stop", "reference", "raw_result"
    ]
    cols = [c for c in preferred if c in df.columns] + [
        c for c in df.columns if c not in preferred
    ]
    return df[cols]


def resistance_explorer(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Search across samples by antibiotic/group, class, subclass or determinant."""
    target = str(target or "").strip().lower()
    if not target:
        return df.iloc[0:0].copy()

    mask = pd.Series(False, index=df.index)
    for col in [
        "Associated antibiotic/group", "class", "subclass",
        "determinant", "element_name", "organism", "sample"
    ]:
        if col in df.columns:
            mask |= df[col].astype(str).str.lower().str.contains(target, regex=False, na=False)

    return df.loc[mask].copy()


def sample_profile(df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """Return all genomic findings belonging to one sample."""
    if "sample" not in df.columns:
        return df.iloc[0:0].copy()
    return df[df["sample"].astype(str).str.lower() == str(sample_id).strip().lower()].copy()
