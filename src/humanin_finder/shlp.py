"""SHLP annotation module for Mitochondrial-Derived Peptides (SHLP1–6).

Based on Cobb et al. 2016 (Aging, PMID: 27070352). SHLP loci are located by their position in
MT-RNR2 (multi-clade 16S locus profile, ``core.locate_mdp_loci``), as Humanin is, and each frame is
read as the species' own ORF in the peptide's primary translation context (Table 1 cytosolic or
Table 2 mitochondrial). Evidence tiers follow ``shlp_references.json`` (E0–E4).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .core import STATUS_COMPLETE, locate_mdp_loci, mdp_display_name

DATA_DIR = Path(__file__).resolve().parent / "data"
SHLP_REF_PATH = DATA_DIR / "shlp_references.json"


def load_shlp_references() -> dict[str, Any]:
    """Load reference sequences and coordinates for SHLP1–6."""
    if not SHLP_REF_PATH.exists():
        return {}
    with open(SHLP_REF_PATH) as f:
        return json.load(f)


def annotate_shlps(
    target_dna: str,
    locus_start_offset: int = 0,
    min_score: float | None = None,
    profiles_path: str | Path | None = None,
    circular: bool = False,
) -> list[dict[str, Any]]:
    """Locate SHLP1–6 by their MT-RNR2 position in a sequence (16S locus, mtDNA, contig or chromosome).

    Every 16S copy at or above the profile's evidence threshold yields one frame per SHLP, with its
    structural status (``Complete``/``Incomplete``/``Interrupted``) and its ungapped identity to the
    human peptide (descriptive only: SHLP loci lie inside the 16S rRNA gene, so peptide similarity
    follows rRNA conservation). The literature tier of the human peptide (``reference_evidence_tier``)
    is carried over only for a ``Complete`` frame identical to it. Loci that are not vertebrate-like
    (``locate_mdp_loci``) are named ``SHLP<n>_position``: they mark the homologous position only.
    With ``circular`` a locus crossing the origin of the record is read through it (``circular_wrap``).
    """
    shlp_dict = load_shlp_references().get("shlp_peptides", {})
    if not shlp_dict:
        return []
    annotations = []
    for locus in locate_mdp_loci(
        target_dna, peptides=list(shlp_dict), min_score=min_score, profiles_path=profiles_path, circular=circular
    ):
        meta = shlp_dict[locus["name"]]
        identical = locus["status"] == STATUS_COMPLETE and locus["seq"] == meta["peptide_sequence_canonical"]
        annotations.append(
            {
                "name": mdp_display_name(locus["name"], locus["vertebrate_like"]),
                "status": locus["status"],
                "seq": locus["seq"],
                "length": locus["length"],
                "strand": locus["strand"],
                "frame": locus["frame"],
                "start": locus["start"] + locus_start_offset,
                "end": locus["end"] + locus_start_offset,
                "circular_wrap": locus["circular_wrap"],
                "translation_table": locus["table"],
                "synteny_score": locus["synteny_score"],
                "start_offset_codons": locus["start_offset_codons"],
                "locus_profile": locus["locus_profile"],
                "locus_density": locus["locus_density"],
                "vertebrate_like": locus["vertebrate_like"],
                # In silico homology evidence for this locus (E1); the literature tier describes the
                # human reference peptide and is only carried over for identical complete frames
                "evidence_tier": "E1",
                "reference_evidence_tier": meta.get("cobb_2016_evidence_level") if identical else None,
                "identity_to_reference": locus["identity_to_reference"],
            }
        )
    return sorted(annotations, key=lambda a: (a["start"], a["name"]))
