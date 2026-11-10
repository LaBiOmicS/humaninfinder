# HumaninFinder

Homology-based location and annotation of **Humanin** and **SHLP1-6**, the mitochondrial-derived
peptides encoded inside the 16S rRNA gene (`MT-RNR2`), in mitogenomes, contigs and nuclear assemblies
(NUMTs).

HumaninFinder finds each peptide by its **position in MT-RNR2**, using multi-clade vertebrate profile
HMMs of the gene. It then reads the species' own reading frame at that position. Every reported frame
states which evidence supports it. All evidence is in silico homology; HumaninFinder makes no claim
that a peptide is expressed or functional.

> Status: **1.0.0b1 (beta)**. It accompanies a manuscript in preparation; see [Citation](#citation).

## Installation

```sh
pip install humanin-finder            # from PyPI
# or, from a clone of this repository
pip install .
```

Requires Python ≥ 3.10, Biopython, click and pyhmmer (HMMER runs in memory; no external binaries are
needed). Biopython and pyhmmer are pinned below their next major version on purpose: the reference
conservation figures depend on them (see `tests/test_conservation_regression.py`).

Check that the installation and every bundled data file are usable:

```sh
humanin-finder setup
```

## Quick start

```sh
# Human reference mitogenome (bundled), vertebrate mitochondrial code
humanin-finder predict -i examples/rcrs_nc_012920_1.fasta -o rcrs -g 2 --circular
```

This writes `rcrs_results.{csv,fasta,gff3,bed}` and `rcrs_run.json`. It reports seven frames:
Humanin (`MAPRGFSCLLLLTSEIDLPVK`, 21 aa under the mitochondrial code) and SHLP1-6, all `Complete`.

## Commands

| Command | Use |
|---|---|
| `predict` | Best Humanin frame per genetic code, plus SHLP1-6, for each record. The usual choice for mitogenomes. |
| `scan` | Every non-redundant locus per record: mitochondrial and NUMT copies in contigs or chromosomes. |
| `conservation` | Coding-evidence controls over a cohort: does conservation at the Humanin window follow the peptide or the host 16S rRNA? |
| `lint` | Validate a FASTA file without changing it. |
| `setup` | Verify the installation and the bundled data files. |

`predict` and `scan` share these options:

- `-i/--input`: FASTA, plain or gzip-compressed. Run through the linter first; an invalid file stops the run.
- `-o/--output`: output prefix.
- `-g/--table`: NCBI genetic code(s), e.g. `2` or `2,1`. The default is `auto`: Table 2 for mitochondrial
  records (the header says *mitochondri*, or the record is ≤ 50 kb), Table 1 otherwise. `auto` only
  covers vertebrates; other mitogenomes need an explicit code (e.g. `5` for insects, `9` for echinoderms).
- `--circular`: treat records ≤ 50 kb as circular molecules. A locus crossing the origin is then read
  through it. Use it for complete mitogenomes, particularly de novo assemblies, whose origin is arbitrary.
  Longer records are always linear.
- `--no-hmm`, `--no-rescue`, `--no-shlp`, `--no-synteny`: turn off the PF15040 domain search, the
  HMM-guided rescue of interrupted frames, the SHLP annotation, or the MT-RNR2 locus search.
- `--all-candidates` (`predict`) / `--best-only` (`scan`): report every candidate, or only the best one per code.

## How a frame is called

1. **Locus.** The record is searched with `nhmmer` against two profiles: the whole MT-RNR2 gene and
   the Humanin core. Each profile hit at or above the profile's evidence threshold maps the start codon
   of every peptide. The threshold is calibrated above the paralogous bacterial and plastid 23S and
   nuclear 28S rRNAs.
2. **Frame.** The species' own ORF is read from the start codon nearest the mapped one (within 6
   codons either way). Humanin is read under the requested code. SHLP1-6 are read under the
   cytosolic code (Table 1) in every record: SHLP1-4 give their published sequences only under it,
   and SHLP5-6, identical under Tables 1 and 2 in human, follow the family. Cobb et al. 2016 do not
   name a code; the reasoning is recorded in `shlp_references.json` (`translation_table_note`).
3. **Status** (structure only, no functional claim):
   - `Complete`: start codon through an in-frame stop. The start is the species' own and may sit
     up to 6 codons from the human one (`start_offset_codons`), with any initiator of the code, so
     `Complete` does not mean the human start codon is conserved.
   - `Incomplete`: no start codon, or no stop within 60 codons.
   - `Interrupted`: an in-frame stop inside the homologous span under the evaluated code.
4. **Name** (how far the evidence goes):
   - `Humanin` / `SHLP<n>`: at a vertebrate-like MT-RNR2 locus. For Humanin, the name is also given to a
     locus whose frame carries the PF15040 domain at the Pfam gathering threshold.
   - `Humanin_position` / `SHLP<n>_position`: at an MT-RNR2 locus outside the vertebrate pattern. Only
     the homologous position is marked.
   - `Humanin_domain`: a PF15040 domain outside any MT-RNR2 locus. Only the domain is supported.

With several codes (`-g 2,1`) names carry a suffix such as `_T2`.

A peptide can be absent from the table even though its locus is present: its start codon may fall in a
deletion or an assembly gap, leaving no frame to read. Those cases are listed in `<prefix>_gaps.csv`
and announced at the end of the run.

## Outputs

| File | Content |
|---|---|
| `<prefix>_results.csv` | One row per reported frame (columns below). |
| `<prefix>_results.fasta` | Peptides; header `sample|name|status|coords_<start>-<end>|evidence_<source>`. |
| `<prefix>_results.gff3` | `ORF` features, 1-based inclusive. Circular records get a `region` with `Is_circular=true`. |
| `<prefix>_results.bed` | 0-based half-open. A feature crossing the origin is split into two blocks with the same name. |
| `<prefix>_gaps.csv` | Peptides whose start codon aligns to a gap (only written when there are any). |
| `<prefix>_run.json` | Provenance: command, options, input SHA-256, version and row counts. |

Main CSV columns:

- **Location:** `start_0based`/`end_0based` (half-open) and `start_1based`/`end_1based` (inclusive)
  include the stop codon of a `Complete` frame. When `circular_wrap` is `True` the feature crosses the
  origin: `end` exceeds the record length, following the GFF3 convention for circular sequences.
- **MT-RNR2 locus:** `synteny_score` is the bitscore of the locus profile. `locus_profile` names that
  profile, and `locus_density` gives bits per aligned profile column.
- **PF15040 domain:** `hmm_bitscore` and `hmm_supported` (at or above the gathering threshold).
  `hmm_evalue` is the E-value of a single comparison (Z = 1), so the same domain gets the same value in
  a short or a long record.
- **Evidence:** `evidence_source` (`Synteny`, `HMM` or both) and `evidence_tier` (`E1`, in silico
  homology). `reference_evidence_tier` is the literature tier of the human peptide (Cobb et al. 2016:
  E2 transcript, E3 protein detected, E4 exogenous activity; SHLP1 E3, SHLP2-4 and SHLP6 E4, SHLP5 E2),
  given only for a `Complete` frame identical to it.
- **Comparison with human:** `identity_to_reference` is the ungapped identity of the homologous span to
  the human peptide. It is descriptive only and does not enter the status.
- **Reading:** `start_offset_codons` is where the local start lies relative to the mapped human start.
  `genetic_code` is the code the frame was read under. `table_1_pep` is the cytosolic reading of a
  Table 2 sORF that starts with a Table 1 initiator.

## Conservation controls

```sh
humanin-finder conservation -i cohort.fasta -o cohort --workers 8              # different species
humanin-finder conservation -i humans.fasta.gz -o humans --population --workers 8   # one species
```

This command projects each genome onto the bundled rCRS genes (MT-RNR2 and MT-CO1). It does **not** use
the locus profiles or the Pfam domain, so its result is independent of the locus caller. For the
Humanin window (rCRS 2633-2695, 21 codons) it reports:

- the share of variation on the third codon position, in the peptide frame and in the two shifted frames;
- amino-acid identity in the same three frames;
- the rank of each statistic among every window of equal length across the 16S gene;
- motif ranks for `CLLL` and `TSEIDLPVK`;
- whether the downstream segment relaxes in genomes with an upstream stop.

MT-CO1 is the positive control. The full record goes to `<prefix>_conservation.json`.

Limits:

- Variation is counted per site, without phylogenetic correction.
- The null windows overlap, so the ranks describe the data and are not a test.
- A 21-codon window has little power. A negative result means "no signal distinguishable from the rRNA
  background", not "no selection".

## Limitations

- Evidence is homology only: position in MT-RNR2 and/or the PF15040 domain. Expression and function
  have to be shown experimentally.
- The profiles were built from vertebrate mitogenomes. Outside vertebrates a locus is named `_position`
  at most.
- `auto` treats records ≤ 50 kb as mitochondrial. A NUMT inside a short nuclear contig needs `-g 1`.
- The PF15040 domain search itself is linear. A Humanin locus crossing the origin still gets its
  PF15040 score, which is computed from the peptide that was read.

## Development

```sh
pip install -e . pytest ruff
ruff check src tests
pytest                     # a few minutes; most of it is the conservation tests
```

To run the regression test on the 343-primate cohort, set `HUMANINFINDER_PRIMATES343` to its FASTA file.

## Citation

If you use HumaninFinder, cite the paper and the version you ran; see [`CITATION.cff`](CITATION.cff).
The manuscript is in preparation.

## License

MIT; see [`LICENSE`](LICENSE). The bundled PF15040 profile comes from Pfam (CC0). The rCRS reference
(NC_012920.1) is from NCBI RefSeq.
