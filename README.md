> [!WARNING]
> **The relaxed-selection results below are withdrawn, and sequence from use-restricted
> genomes has been removed from this repository. See [NOTICE.md](NOTICE.md) (2026-08-28).**
> The gene-count results are unaffected.

# Relaxed selection and convergent loss of PCWDE in ectomycorrhizal fungi

Analysis code, intermediate data, results, and figures for:

> **Relaxed selection and convergent loss erode the retained plant-cell-wall-degrading enzymes of ectomycorrhizal fungi.**
> Minseo Kim, Jae-Ho Shin. *Preprint / in preparation.* DOI: _[TBD — bioRxiv]_

This repository reproduces the study's two central analyses on 182–183 fungal genomes:
- **H1 (gene-level selection):** HyPhy **RELAX** on PCWDE codon alignments, ECM lineages as test branches.
- **H3 (species-level gene loss):** **PGLS** of PCWDE gene counts vs. lifestyle, phylogeny-corrected.

~~Main finding: gene loss and relaxed selection are **separable axes**. The dedicated cellulose/pectin
families (GH6, GH7, GH28) are relaxed on ECM branches (K = 0.46–0.87) *and* largely lost, whereas the
multifunctional GH5 is reduced in copy number but kept under purifying selection (K = 1.00).~~

**Withdrawn (2026-08-28).** The relaxed-selection half of this statement does not survive
re-analysis; see [NOTICE.md](NOTICE.md). The gene-loss half stands.

## Repository layout

```
scripts/     analysis pipeline (Python, R, shell)
data/        inputs & intermediates
  genome_accessions.tsv   183 genomes: JGI MycoCosm portal, species, lifestyle, order
  pcwde_counts.tsv        per-genome PCWDE / total CAZyme counts (dbCAN3)
  pgls_data.tsv           per-species table for PGLS (counts + functional counts)
  perfamily_intact.tsv    per-species intact counts per family (GH5/6/7/28)
  species_tree.nwk        OrthoFinder species tree (182 tips)
  calls/                  catalytic-residue calls per family (intact/disrupted/fragment)
  alignments/             RELAX codon alignments, published genomes only (see NOTICE.md) + labelled trees
results/
  RELAX_RESULTS.md        H1 write-up (K, P, LRT per family)
  PGLS_RESULTS.md         H3 write-up (effect sizes, Pagel's λ, per-family bridge)
  relax_json/             raw HyPhy RELAX output per family
  pgls_output.txt, perfamily_output.txt   raw R output
figures/     Figure 1–4 (PDF vector + PNG)
```

Raw genomes/proteomes are **not** redistributed here; they are public at JGI MycoCosm and NCBI
(see `data/genome_accessions.tsv`).

## Pipeline overview

**Phase 1 — assembly & annotation.** Genomes/proteomes from JGI MycoCosm + NCBI → header cleanup
(`clean_proteomes.py`) → BUSCO QC (fungi_odb10 ≥ 85 %) → OrthoFinder species tree → dbCAN3 CAZyme
annotation (`annotate_dbcan.sh`) → PCWDE counts (`extract_pcwde.py`).

**Phase 2 — H1 (selection).**
1. Catalytic-residue scoring vs. characterized references (`catalytic_check.py`, `run_catalytic_all.sh`).
2. CDS retrieval (`download_cds.py`) and codon alignment (`codon_align.py`; maps proteinId↔transcriptId
   via the shared JGI model name), column trimming (`trim_codon.py`), GH5 subsampling
   (`subsample_codon.py`).
3. Gene tree (FastTree) → label ECM tips (`label_tree_trim.py`) → HyPhy RELAX (`run_relax_trim.sh`,
   `run_relax_sub.sh`). Parse with `parse_relax.py`, summarize with `summarize_relax.py`.

**Phase 3 — H3 (PGLS).** Build tables (`build_pgls_table.py`, `build_perfamily_table.py`) → PGLS and
phylogenetic logistic regression in R (`pgls_analysis.R`, `pgls_perfamily.R`).

**Figures.** `figures.R`.

## Reproducing

Software: Python 3.10, R 4.5 (`ape`, `caper`, `phylolm`, `nlme`, `geiger`, `ggplot2`, `ggtree`,
`patchwork`, `ggrepel`), HyPhy 2.5, FastTree, MAFFT, OrthoFinder, BUSCO, run_dbcan (dbCAN3).

From the intermediate data provided here, the two headline analyses can be rerun directly:

```bash
# H1 — RELAX (per family)
for fam in GH6 GH7 GH28 GH5; do bash scripts/run_relax_trim.sh $fam 0.5; done
python3 scripts/summarize_relax.py

# H3 — PGLS  (conda env with r-ape/caper/phylolm)
Rscript scripts/pgls_analysis.R
Rscript scripts/pgls_perfamily.R

# Figures
Rscript scripts/figures.R
```

Paths inside the scripts point to the original project directory (`~/projects/aicomp-h1h3`); adjust the
`PROJ` variable at the top of each script to this repository root to rerun in place.

## Key results

| Family | Enzyme | RELAX K (H1) | ECM retention (H3) |
|---|---|:--:|---|
| GH7 | cellobiohydrolase I | 0.46 (P = 1.6e-9) | present in 13 % of ECM |
| GH6 | cellobiohydrolase II | 0.80 (P = 1.8e-4) | present in 11 % of ECM |
| GH28 | polygalacturonase | 0.87 (P = 2.7e-5) | 2.1 vs 6.4 copies |
| GH5 | multifunctional GH | 1.00 (n.s.) | 100 % of ECM, 12.0 vs 19.4 copies |

PGLS: ECM carry ≈0.41× the raw and ≈0.55× the functional PCWDE of free-living decomposers
(P ≤ 9.8e-10; Pagel's λ ≈ 0.93–1.0), i.e. a convergent reduction across 11 fungal orders.

## Citation

If you use this code or data, please cite the paper (above) and this repository:

> Kim M. (2026). Analysis code and data for "Relaxed selection and convergent loss erode the retained
> PCWDE of ectomycorrhizal fungi". Zenodo. DOI: _[TBD]_

## Contact

Minseo Kim — 2023024947@knu.ac.kr · Department of Applied Biosciences, Kyungpook National University.

## License

Code: MIT (see `LICENSE`). Data, results, and figures: CC-BY-4.0.
