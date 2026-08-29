# NOTICE — withdrawn claims and redacted data (2026-08-28)

Two corrections apply to this repository. Both are described here rather than made silently.

## 1. The relaxed-selection claims of the v1.0 release are withdrawn

The v1.0 release and the accompanying preprint reported significant relaxed selection on
ectomycorrhizal branches for *GH6*, *GH7* and *GH28*. **Those claims are withdrawn.** On
re-analysis the branch test used to obtain them does not hold its nominal error rate under
the labelling design applied here, and the fits are not reproducible across optimiser
configurations. The gene-count result — ectomycorrhizal genomes carry markedly fewer
plant-cell-wall-degrading enzyme genes than free-living decomposers — is unaffected and
has been reproduced.

The successor manuscript sets out the grounds in full. Until it appears, the safe reading
of this repository is that its count-based results stand and its selection-based results
do not.

## 2. Sequence from use-restricted genomes has been removed

The v1.0 release included coding sequence from 16 JGI MycoCosm genomes that have no
associated publication and may be use-restricted. We had not established their
use-restriction status before that release was published. Removed on 2026-08-28:

| What | Where | Removed |
|---|---|---|
| Codon alignments | `data/alignments/*.codon.fasta` | 437 sequences from the 16 genomes |
| Per-taxon codon states | `results/relax_json/*.json`, `substitutions` block | all |
| Per-site log-likelihoods | `results/relax_json/*.json`, `Site Log Likelihood` block | all |

Alignments now contain only genomes with an associated publication. Model fits, test
results, branch attributes, trees, gene counts, catalytic-residue calls, scripts and
figures are unchanged: every reported statistic in the v1.0 release can still be checked
against the retained files, and the alignments themselves can be rebuilt from the JGI
MycoCosm portal by anyone entitled to the underlying data.

**What has and has not been removed.** We would rather set this out than let a reader assume
more than is true.

- **GitHub — done.** Rewriting history was not sufficient: the earlier objects continued to be
  served when addressed by their commit SHA, and that SHA was not secret, since it appears in
  the directory name inside the v1.0 Zenodo archive. The repository was therefore deleted and
  recreated from the redacted tree. The earlier commit, its files and the v1.0 tarball now
  return HTTP 404.
- **Zenodo — partly.** A redacted version is deposited as doi:10.5281/zenodo.22146740. The
  original deposit, doi:10.5281/zenodo.21449619 of 20 July 2026, **is still downloadable**: a
  Zenodo record cannot be withdrawn by its depositor, and we have asked Zenodo Support to
  restrict its files while keeping the record and DOI in place so the correction stays visible.
  Until they act, the sequence remains obtainable there.
- **Already downloaded copies are beyond our reach.** The v1.0 archive was public under CC-BY
  for five weeks and recorded two downloads.

We are writing to the JGI Data Portal to report this and to ask which of the genomes are
in fact use-restricted.

Questions: MinSeo Kim <2023024947@knu.ac.kr>, Prof. Jae-Ho Shin <jhshin@knu.ac.kr>.
