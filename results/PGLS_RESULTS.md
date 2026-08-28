# H3 PGLS — PCWDE reduction tracks ECM lifestyle after phylogenetic correction

**Date:** 2026-07-17 · Phase 3 (H3) — PGLS complete.

## Question (H3)
Does the reduced plant-cell-wall-degrading-enzyme (PCWDE) repertoire of ECM fungi track lifestyle
**after** correcting for phylogeny? Miyauchi (2020) showed ECM have fewer PCWDE (non-phylogenetic
Mann-Whitney, here re-confirmed p=1.1e-12), but closely related species share counts by descent, so the
rigorous test is a **phylogenetic generalized least squares (PGLS)** regression.

## Data
- **Tree**: OrthoFinder species tree (`SpeciesTree_rooted.txt`, 182 tips = QC-passed genomes), node labels
  dropped, polytomies resolved (`multi2di`).
- **Response(s)**: `log10(count+1)` of
  1. **raw PCWDE count** — all dbCAN PCWDE families (`annotation/pcwde_counts.tsv`).
  2. **functional (intact) count** — catalytically-intact copies in GH5/6/7/28 (from Step-1 `*_calls.tsv`;
     AA9 excluded, SignalP-dependent). This is the novel "capacity" response the count alone misses.
- **Predictor**: `is_ECM` (ECM vs other).
- **Method**: `caper::pgls` with Pagel's λ estimated by ML; ordinary `lm` (OLS) reported alongside as the
  no-phylogeny baseline. Table builder `scripts/build_pgls_table.py` → `h3_pgls/pgls_data.tsv`.

## Results

| contrast | response | ECM effect (log10) | ECM ×count | **PGLS p** | Pagel λ | OLS p (no phylo) |
|---|---|:--:|:--:|:--:|:--:|:--:|
| **A) ECM vs SAP+WD** (n=122 vs 33) | raw PCWDE | −0.387 | **0.41×** | **1.6e-11** | 0.94 | 1.6e-19 |
| A) ECM vs SAP+WD | functional intact | −0.262 | **0.55×** | **9.8e-10** | 0.93 | 5.8e-15 |
| **B) ECM vs all non-ECM** (n=122 vs 60) | raw PCWDE | −0.362 | **0.43×** | **2.0e-8** | 1.00 | 1.3e-19 |
| B) ECM vs all non-ECM | functional intact | −0.237 | **0.58×** | **7.1e-7** | 1.00 | 3.5e-16 |

(SAP = soil saprotroph, WD = wood decayer — the free-living decomposers, i.e. the Miyauchi contrast.)

## Interpretation (H3 answer)
- ✅ **H3 confirmed**: ECM lifestyle predicts a large, highly significant PCWDE reduction **even after
  phylogenetic correction** (all p ≤ 7e-7). ECM species carry ~0.4× the raw PCWDE and ~0.55× the functional
  count of free-living decomposers.
- **Not phylogenetic inertia**: Pagel's λ ≈ 0.93–1.0 means counts carry strong phylogenetic signal (relatives
  resemble each other), yet the ECM effect survives. Because ECM is polyphyletic (122 species across
  Boletales/Russulales/Agaricales/Pezizales — ~12 independent origins), PGLS treats those origins as
  replication: the reduction is **convergent**, repeated at each independent ECM origin, not a one-clade
  artifact. PGLS shrinks the effect vs OLS (e.g. −0.39 vs −0.45 log10) but keeps it significant — the honest,
  phylogenetically-controlled version of Miyauchi's pattern.
- **Ties to H1**: H1 (RELAX) showed the PCWDE that ECM *retain* are under relaxed selection for the dedicated
  cellulose/pectin families (GH6/7/28), maintained only for multifunctional GH5. H3 (PGLS) shows ECM also
  *shed* PCWDE genes wholesale. Together: **ECM both lose PCWDE genes (H3, count) and relax selection on the
  dedicated ones they keep (H1, selection)** — a two-front erosion of lignocellulolytic capacity that the raw
  gene count only partially captures.

## Per-family bridge (H1 ↔ H3) — the payoff

Testing each catalytic family separately (ECM vs SAP+WD, phylogeny-corrected) against its RELAX K from H1.
GH7/GH6 are near presence/absence in ECM → phylogenetic logistic regression (`phyloglm`); GH28/GH5 are
prevalent → log-count PGLS.

| family | RELAX K (H1) | H3 test | ECM effect | **p (phylo)** | ECM vs DEC |
|---|:--:|---|:--:|:--:|---|
| **GH7** | 0.46 (most relaxed) | phyloglm presence | −1.76 logit | 8.7e-4 | present in **13% vs 76%** |
| **GH6** | 0.80 | phyloglm presence | −2.53 logit | 1.3e-4 | present in **11% vs 79%** |
| **GH28** | 0.87 | log-count PGLS | −0.251 | 2.9e-5 | 2.1 vs 6.4 copies |
| **GH5** | 1.00 (maintained) | log-count PGLS | −0.192 | 1.3e-7 | 12.0 vs 19.4 copies |

**Two separable axes.** Count reduction (H3) hits *every* family in ECM (all p significant) — even the
selection-maintained GH5 is count-reduced (12 vs 19). But selection relaxation (H1) hits only the *dedicated*
families. So:
- **GH7/GH6** (dedicated cellulose CBH): **lost outright** (present in ~12% of ECM) **and** relaxed (K<1) —
  double erosion.
- **GH5** (multifunctional): **count-reduced but selection-maintained** (K=1.0) — the copies ECM keep are
  still doing real work (likely non-PCWDE roles).
- The RELAX-K → retention gradient is monotone-ish: K 0.46/0.80 → ~12% kept; K 0.87 → 74%; K 1.0 → 100%.

→ "Reduced count" and "relaxed selection" are **not the same thing**, and the raw gene count conflates them.
The dedicated lignocellulose machinery is being erased on both axes; the multifunctional family survives
because its other jobs keep it under selection despite fewer copies.

## Caveats
- Functional count covers 4 assessed families (GH5/6/7/28), not the full PCWDE set; a proxy, not the whole
  capacity. AA9 pending SignalP.
- λ estimates at the boundary (1.0) in contrast B → ML at the Brownian limit; effect is robust regardless.
- log10(count+1) with a binary predictor (standard for gene-family sizes); count-specific models (phylo
  Poisson) not used.
- Response is count/functional-count; "conditional saprotroph" ecology/expression evidence (the within-ECM
  reading of H3) is not tested here (no expression data).

## Next / possible extension
- Phylogenetic ANOVA across the full lifestyle gradient (rare categories permitting).
- Within-ECM "conditional saprotroph" test needs expression/ecology data (not in hand).
- AA9 family added once SignalP resolves its secretome numbering.

## Files (lab101 `~/projects/aicomp-h1h3/h3_pgls/`)
- `pgls_data.tsv`, `perfamily_intact.tsv` (analysis tables); `pgls_output.txt`, `perfamily_output.txt` (raw model output); this doc.
- Scripts: `scripts/{build_pgls_table.py, build_perfamily_table.py, pgls_analysis.R, pgls_perfamily.R}`. Env: conda `pgls` (r-ape/caper/phylolm).
