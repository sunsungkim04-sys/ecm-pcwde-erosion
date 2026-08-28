#!/usr/bin/env Rscript
# H1<->H3 bridge: per-family phylo-corrected test of ECM reduction, vs the RELAX K from H1.
# GH7/GH6 are near presence/absence in ECM -> phylogenetic logistic regression (phyloglm).
# GH28/GH5 are prevalent -> log-count PGLS. All vs SAP+WD (free-living decomposers).
suppressMessages({library(ape); library(caper); library(phylolm)})
PROJ <- "/home1/minseo1101/projects/aicomp-h1h3"
tr <- read.tree(file.path(PROJ,"orthofinder/run/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt"))
tr$node.label <- NULL; tr <- multi2di(tr); tr$edge.length[tr$edge.length<=0] <- 1e-8
d  <- read.delim(file.path(PROJ,"h3_pgls/perfamily_intact.tsv"), stringsAsFactors=FALSE)
d  <- d[d$lifestyle %in% c("ECM","SAP","WD"), ]           # ECM vs free-living decomposers
d$is_ECM <- ifelse(d$lifestyle=="ECM", 1L, 0L)            # 1=ECM
relaxK <- c(GH7=0.46, GH6=0.80, GH28=0.87, GH5=1.00)

cat("family  RELAX_K  test                 ECM_effect     p          note\n")
for (fam in c("GH7","GH6","GH28","GH5")) {
  dd <- d[d$portal %in% tr$tip.label, ]
  tt <- keep.tip(tr, dd$portal); dd <- dd[match(tt$tip.label, dd$portal),]
  rownames(dd) <- dd$portal
  cnt <- dd[[fam]]
  presrate_ecm <- mean(cnt[dd$is_ECM==1]>0); presrate_dec <- mean(cnt[dd$is_ECM==0]>0)
  if (mean(cnt>0) < 0.5) {
    # presence/absence -> phylogenetic logistic regression
    dd$pres <- as.integer(cnt>0)
    m <- tryCatch(phyloglm(pres ~ is_ECM, data=dd, phy=tt, method="logistic_MPLE", btol=30),
                  error=function(e) NULL)
    if (is.null(m)) { cat(sprintf("%-6s  %5.2f   phyloglm(pres)       FAILED (separation); ECM %.0f%% vs DEC %.0f%% present\n",
                                  fam, relaxK[fam], 100*presrate_ecm, 100*presrate_dec)); next }
    s <- summary(m); co <- s$coefficients
    cat(sprintf("%-6s  %5.2f   phyloglm(presence)   %+.2f (logit)   p=%.3g   ECM %.0f%% vs DEC %.0f%% present\n",
                fam, relaxK[fam], co["is_ECM","Estimate"], co["is_ECM","p.value"], 100*presrate_ecm, 100*presrate_dec))
  } else {
    # count -> log-count PGLS (Pagel lambda ML)
    dd$y <- log10(cnt+1)
    cd <- comparative.data(tt, dd, names.col="portal", vcv=TRUE, warn.dropped=FALSE)
    m  <- pgls(y ~ is_ECM, cd, lambda="ML"); s <- summary(m); co <- s$coefficients
    lam <- tryCatch(s$param.CI$lambda$opt, error=function(e) NA)
    cat(sprintf("%-6s  %5.2f   log-count PGLS       %+.3f (log10)  p=%.3g   ECM %.1f vs DEC %.1f mean; lambda=%.2f\n",
                fam, relaxK[fam], co["is_ECM","Estimate"], co["is_ECM","Pr(>|t|)"],
                mean(cnt[dd$is_ECM==1]), mean(cnt[dd$is_ECM==0]), lam))
  }
}
cat("\nGradient: strongest RELAX relaxation (GH7 K=0.46) = most reduced/lost in ECM; maintained (GH5 K=1.0) = kept.\n")
