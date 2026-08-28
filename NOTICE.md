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
associated publication and may be use-restricted. It should not have been redistributed
before their portal status was checked. Removed on 2026-08-28:

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

**This removal is not yet complete, and we would rather say so than imply otherwise.**

- The v1.0 tag, its release and its source tarball have been deleted, and the repository
  history was rewritten, so nothing here now links to the withdrawn files.
- **The pre-rewrite objects are nonetheless still served by GitHub when addressed by their
  commit SHA**, because unreferenced objects are not garbage-collected on request. We have
  asked GitHub Support to purge them; until that is done, the earlier files remain
  reachable by anyone who knows the SHA, and the SHA is not secret — it appears in the
  directory name inside the Zenodo archive.
- The v1.0 release was archived on Zenodo (doi:10.5281/zenodo.21449619) on 20 July 2026
  under CC-BY, and a Zenodo record cannot be withdrawn. A redacted version has been
  deposited and Zenodo Support asked to restrict the files of the earlier one.
- Anything already downloaded is beyond our reach.

We are writing to the JGI Data Portal to report this and to ask which of the genomes are
in fact use-restricted.

Questions: MinSeo Kim <2023024947@knu.ac.kr>, Prof. Jae-Ho Shin <jhshin@knu.ac.kr>.
