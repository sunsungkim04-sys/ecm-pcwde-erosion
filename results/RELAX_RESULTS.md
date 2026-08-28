# H1 RELAX — Relaxed selection on retained PCWDE in ECM fungi

**Date:** 2026-07-16 · Phase 2 Step 2 (RELAX) complete.

## Question
Are the plant-cell-wall-degrading enzyme (PCWDE) genes that ECM fungi *retain* under maintained purifying
selection (→ functional, "conditional saprotroph") or under **relaxed** selection (→ eroding relics)?
Step 1 showed their catalytic residues are 86–100 % intact — but intactness ≠ maintained selection.
HyPhy **RELAX** tests this directly: ECM tips = **Test**, all other lifestyles = **Reference**; K<1 & p<0.05 = relaxed.

## Pipeline (what was run)
1. **CDS restore + download** — JGI tape-restored GeneCatalog CDS; `download_cds.py` → **51/52** species on disk
   (1 ECM sp. *Gaumor1_1* still `RESTORE_REGISTERED`, deferred). All gzip-valid.
2. **CDS↔protein ID-mapping fix (critical)** — JGI labels the proteome by *proteinId* but the CDS by *transcriptId*;
   the two numeric IDs differ per release (only the shared model name, 4th `|` field, is stable). The old
   `codon_align.py` keyed on the numeric field → silently dropped 50–85 % of genes for version-mismatched species
   (Phchr2 13 %, Thacu1 23 % mapped) and kept ~100 % only where IDs coincided (Lacbi2). Rewritten to bridge
   proteinId → model name → nt (fallback: translated-sequence match). Recovery → **100 %** for all species tested.
3. **codon_align recovery** (intact proteins back-translated to codon alignments):

   | family | before fix | after (seqs) | ECM tips | len_mismatch |
   |---|---|---|---|---|
   | GH7 | 4 (ECM 0) | 140 | 13 | 0 |
   | GH6 | — | 69 | 12 | 0 |
   | GH28 | 15 (ECM 3) | 289 | 44 | 0 |
   | GH5 | — | 902 | 255 | 2 |
4. **Column trimming** — alignments were gap-dominated (GH5 18 943 codons, only ~350 informative). Dropped
   codon columns with >50 % gaps (`trim_codon.py`, maxgap 0.5) → 353–507 informative codons/family. Standard
   practice; also required for tractability (untrimmed GH7 did not finish in >1 h).
5. **GH5 subsample** — 902 tips is intractable for RELAX (branch-dominated). Subsampled to 80 ECM + 80 reference
   (`subsample_codon.py`, seed 42). Full 902-tip GH5 left running as confirmation.
6. **RELAX** — per family: FastTree (`-nt -gtr`) → label ECM `{Test}` → `hyphy relax`. K & p from JSON `test results`.

## Results

| family | enzyme (function) | nseq | ECM(Test) | sites | **K** | **p** | LRT | verdict |
|---|---|---|---|---|---|---|---|---|
| **GH7** | cellobiohydrolase I (crystalline cellulose) | 140 | 13 | 507 | **0.46** | 1.6e-9 | 36.4 | **relaxed** |
| **GH6** | cellobiohydrolase II (cellulose) | 69 | 12 | 397 | **0.80** | 1.8e-4 | 14.0 | **relaxed** |
| **GH28** | polygalacturonase (pectin) | 289 | 44 | 376 | **0.87** | 2.7e-5 | 17.6 | **relaxed** |
| **GH5** (subsample) | broad GH (cellulase/mannanase/β-glucanase; incl. fungal cell wall) | 160 | 80 | 353 | 0.96 | 0.14 | 2.2 | maintained (ns) |
| **GH5** (full) | " | 902 | 255 | 353 | **1.00** | **1.0** | −0.07 | **maintained (ns)** |

GH5 full (all 255 ECM, 902 tips) **confirms** the subsample: K=1.00, p=1.0, and the Reference/Test ω
distributions are identical — no relaxation whatsoever.

ω distributions (RELAX alternative, Reference → Test): the high-ω tail is damped on ECM branches for the relaxed
families — GH7 19.5→3.97, GH28 11.0→8.0 — while GH5 is essentially unchanged (41.7→36.0).

## Interpretation (H1 answer)
**The PCWDE genes ECM fungi retain are NOT uniformly maintained.** The three *dedicated* plant-cell-wall families —
crystalline-cellulose cellobiohydrolases (GH6, GH7) and pectin polygalacturonase (GH28) — are under **significantly
relaxed selection** on ECM branches, i.e. functional erosion is in progress, *despite intact catalytic residues*
(Step 1). Only the **broad, multifunctional GH5** (much of which acts on fungal cell wall / β-glucans, not just
plant cell wall) remains under purifying selection.

→ Directly supports the mother-plan thesis **"Gene Count ≠ Functional Capacity"**: raw counts and even
catalytic-residue intactness overstate the retained lignocellulolytic capacity of ECM fungi; selection signatures
reveal the canonical decay machinery is decaying.

## Caveats
- Family-level tests; GH5 lumps many subfamilies, so "maintained overall" does not exclude relaxation in specific
  subfamilies. GH7/GH6 relaxation is the cleanest (dedicated function).
- GH5 "maintained" is confirmed by the full 902-tip run (K=1.00, p=1.0), not only the subsample.
- Column trimming (maxgap 0.5) applied uniformly; robustness to threshold not yet swept.
- 1 ECM species (*Gaumor1_1*) awaits JGI restore; adds a few GH28/GH5 genes when available.

## Key files (lab101 `~/projects/aicomp-h1h3`)
- Results JSON: `h1_pcwde/codon/{GH7,GH6,GH28}.trim.relax.json`, `GH5.sub.relax.json`
- Alignments: `h1_pcwde/codon/*.trim.codon.fasta` (+ `GH5.sub.codon.fasta`)
- Scripts: `scripts/{codon_align.py (fixed), trim_codon.py, subsample_codon.py, run_relax_trim.sh, run_relax_sub.sh, parse_relax.py, summarize_relax.py}`
- Table: `python3 scripts/summarize_relax.py`
