#!/usr/bin/env Rscript
# H3 PGLS: does lifestyle (ECM vs free-living decomposer) predict PCWDE gene count after phylogenetic correction?
# Responses: log10(raw PCWDE count) and log10(functional=intact count, GH5/6/7/28). Predictor: ECM vs other.
suppressMessages({library(ape); library(caper)})
PROJ <- "/home1/minseo1101/projects/aicomp-h1h3"
tr  <- read.tree(file.path(PROJ,"orthofinder/run/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt"))
dat <- read.delim(file.path(PROJ,"h3_pgls/pgls_data.tsv"), stringsAsFactors=FALSE)

# drop internal node labels (support values collide with tip labels in caper), resolve polytomies, guard zero branches
tr$node.label <- NULL
tr <- multi2di(tr)
tr$edge.length[tr$edge.length <= 0] <- 1e-8

dat$logpcwde <- log10(dat$pcwde_count + 1)
dat$logfunc  <- log10(dat$func_count  + 1)
dat$is_ECM   <- factor(ifelse(dat$lifestyle=="ECM","ECM","other"), levels=c("other","ECM"))  # 'other' baseline

fit <- function(subset_life, resp, label){
  d <- dat[dat$lifestyle %in% subset_life, ]
  d <- d[d$portal %in% tr$tip.label, ]
  tt <- keep.tip(tr, d$portal)
  cd <- comparative.data(tt, d, names.col="portal", vcv=TRUE, warn.dropped=FALSE)
  f  <- as.formula(paste(resp, "~ is_ECM"))
  m  <- pgls(f, cd, lambda="ML")
  s  <- summary(m)
  co <- s$coefficients
  est <- co["is_ECMECM","Estimate"]; p <- co["is_ECMECM","Pr(>|t|)"]
  lam <- tryCatch(s$param.CI$lambda$opt, error=function(e) m$param["lambda"])
  ols <- summary(lm(f, data=d))
  oest<- ols$coefficients["is_ECMECM","Estimate"]; op <- ols$coefficients["is_ECMECM","Pr(>|t|)"]
  nE <- sum(d$is_ECM=="ECM"); nO <- sum(d$is_ECM=="other")
  cat(sprintf("\n=== %s  |  response=%s  (ECM n=%d vs other n=%d) ===\n", label, resp, nE, nO))
  cat(sprintf("  PGLS: ECM effect = %+.3f (log10)  p = %.3g   Pagel lambda = %.3f\n", est, p, lam))
  cat(sprintf("        -> ECM has %.2fx the count of the baseline group (10^effect)\n", 10^est))
  cat(sprintf("  OLS (no phylo): ECM effect = %+.3f  p = %.3g\n", oest, op))
  invisible(list(est=est,p=p,lambda=lam,ols_p=op,nE=nE,nO=nO))
}

cat("################ H3 PGLS ################\n")
cat("Tree tips:", length(tr$tip.label), " Data rows:", nrow(dat), "\n")

# Contrast A: ECM vs free-living decomposers (SAP + WD) — the Miyauchi comparison
fit(c("ECM","SAP","WD"), "logpcwde", "A) ECM vs SAP+WD  [raw PCWDE]")
fit(c("ECM","SAP","WD"), "logfunc",  "A) ECM vs SAP+WD  [functional intact]")

# Contrast B: ECM vs all non-ECM lifestyles
allother <- setdiff(unique(dat$lifestyle), "ECM")
fit(c("ECM", allother), "logpcwde", "B) ECM vs all non-ECM  [raw PCWDE]")
fit(c("ECM", allother), "logfunc",  "B) ECM vs all non-ECM  [functional intact]")

cat("\n################ done ################\n")
