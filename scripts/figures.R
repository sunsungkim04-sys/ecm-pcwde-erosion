#!/usr/bin/env Rscript
# Generate manuscript figures 1-4 for the H1/H3 (ECM PCWDE) paper.
suppressMessages({library(ape); library(ggplot2); library(ggtree); library(patchwork); library(ggrepel); library(scales)})
PROJ <- "/home1/minseo1101/projects/aicomp-h1h3"
OUT  <- file.path(PROJ,"manuscript","figures"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
d  <- read.delim(file.path(PROJ,"h3_pgls/pgls_data.tsv"), stringsAsFactors=FALSE)
pf <- read.delim(file.path(PROJ,"h3_pgls/perfamily_intact.tsv"), stringsAsFactors=FALSE)
tr <- read.tree(file.path(PROJ,"orthofinder/run/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt")); tr$node.label <- NULL

# consistent lifestyle order/colors (ECM emphasized)
life_lv <- c("ECM","SAP","WD","ORM","ERM","ENDO","PATH","AM","YEAST","LICHEN")
pal <- c(ECM="#D7263D", SAP="#1B998B", WD="#2E86AB", ORM="#7B68EE", ERM="#9370DB",
         ENDO="#8FBC8F", PATH="#F46036", AM="#C0C0C0", YEAST="#808080", LICHEN="#A0522D")
d$lifestyle  <- factor(d$lifestyle,  levels=life_lv)
pf$lifestyle <- factor(pf$lifestyle, levels=life_lv)
th <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())

## ---------- Figure 1: species tree + lifestyle + PCWDE count ----------
dd <- d; rownames(dd) <- dd$portal
p1 <- ggtree(tr, size=0.18, layout="rectangular") %<+% dd +
  geom_tippoint(aes(color=lifestyle), size=0.9, na.rm=TRUE) +
  scale_color_manual(values=pal, name="Lifestyle", drop=FALSE) +
  theme_tree2() + theme(legend.position="left")
# aligned PCWDE count bar via facet_plot
bar_df <- data.frame(id=dd$portal, count=dd$pcwde_count)
p1f <- facet_plot(p1, panel="PCWDE count", data=bar_df,
                  geom=geom_segment, aes(x=0, xend=count, y=y, yend=y), linewidth=0.5, color="grey30") +
  theme(strip.background=element_rect(fill="grey92"))
ggsave(file.path(OUT,"Figure1_tree_lifestyle_PCWDE.pdf"), p1f, width=9, height=11)
ggsave(file.path(OUT,"Figure1_tree_lifestyle_PCWDE.png"), p1f, width=9, height=11, dpi=200)

## ---------- Figure 2: RELAX K per family ----------
relax <- data.frame(
  family=c("GH7","GH6","GH28","GH5"),
  enzyme=c("cellobiohydrolase I","cellobiohydrolase II","polygalacturonase","multifunctional GH"),
  K=c(0.46,0.80,0.87,1.00),
  p=c(1.6e-9,1.8e-4,2.7e-5,1.0),
  nECM=c(13,12,44,255))
relax$family <- factor(relax$family, levels=c("GH7","GH6","GH28","GH5"))
relax$sig <- ifelse(relax$p<0.05,"relaxed (K<1, P<0.05)","maintained (n.s.)")
relax$lab <- ifelse(relax$p<0.05, sprintf("K=%.2f\nP=%s", relax$K, format(relax$p, scientific=TRUE, digits=2)),
                    sprintf("K=%.2f\nn.s.", relax$K))
p2 <- ggplot(relax, aes(x=family, y=K, color=sig)) +
  geom_hline(yintercept=1, linetype=2, color="grey50") +
  geom_segment(aes(xend=family, y=1, yend=K), linewidth=0.7) +
  geom_point(aes(size=nECM)) +
  geom_text(aes(label=lab), vjust=1.6, size=3, show.legend=FALSE) +
  scale_color_manual(values=c("relaxed (K<1, P<0.05)"="#D7263D","maintained (n.s.)"="#1B998B"), name=NULL) +
  scale_size_continuous(range=c(3,8), name="ECM test\nbranches") +
  scale_y_continuous(limits=c(0.3,1.12), breaks=seq(0.4,1.0,0.2)) +
  labs(x="PCWDE family (dedicated → multifunctional)", y="RELAX relaxation parameter K",
       title="Selection is relaxed on the dedicated PCWDE of ECM fungi",
       subtitle="K<1 = relaxed selection on ECM branches; dashed line K=1 (no change)") + th +
  theme(legend.position="right")
ggsave(file.path(OUT,"Figure2_RELAX_K.pdf"), p2, width=8, height=5)
ggsave(file.path(OUT,"Figure2_RELAX_K.png"), p2, width=8, height=5, dpi=200)

## ---------- Figure 3: PGLS ----------
# 3a: log PCWDE count by lifestyle (decomposer groups vs ECM)
d3 <- d[d$lifestyle %in% c("ECM","SAP","WD"),]
d3$grp <- factor(ifelse(d3$lifestyle=="ECM","ECM","free-living\ndecomposer (SAP+WD)"), levels=c("free-living\ndecomposer (SAP+WD)","ECM"))
p3a <- ggplot(d3, aes(grp, pcwde_count, fill=grp)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.5) +
  geom_jitter(aes(color=grp), width=0.12, size=1, alpha=0.6, show.legend=FALSE) +
  scale_fill_manual(values=c("ECM"="#D7263D","free-living\ndecomposer (SAP+WD)"="#1B998B")) +
  scale_color_manual(values=c("ECM"="#D7263D","free-living\ndecomposer (SAP+WD)"="#1B998B")) +
  scale_y_log10() + labs(x=NULL, y="PCWDE genes per genome (log)", title="a  Raw PCWDE count") +
  th + theme(legend.position="none")
# 3b: forest plot of ECM effect (fold-change = 10^effect) PGLS vs OLS
eff <- data.frame(
  contrast=rep(c("ECM vs SAP+WD","ECM vs all non-ECM"), each=2),
  response=rep(c("raw PCWDE","functional (intact)"), 2),
  pgls_fold=c(0.41,0.55,0.43,0.58),
  pgls_p  =c(1.6e-11,9.8e-10,2.0e-8,7.1e-7))
eff$lab <- sprintf("%s\n%s", eff$contrast, eff$response)
eff$row <- factor(eff$lab, levels=rev(eff$lab))
p3b <- ggplot(eff, aes(x=pgls_fold, y=row, color=response)) +
  geom_vline(xintercept=1, linetype=2, color="grey50") +
  geom_point(size=3.5) +
  geom_text(aes(label=sprintf("%.2f×  P=%s", pgls_fold, format(pgls_p,scientific=TRUE,digits=2))),
            hjust=-0.15, size=3, show.legend=FALSE) +
  scale_color_manual(values=c("raw PCWDE"="#22223B","functional (intact)"="#D7263D"), name=NULL) +
  scale_x_continuous(limits=c(0.3,1.15), breaks=c(0.4,0.6,0.8,1.0)) +
  labs(x="ECM PCWDE relative to comparator (PGLS, phylogeny-corrected)", y=NULL,
       title="b  ECM reduction survives phylogenetic correction (λ≈0.93–1.0)") +
  th + theme(legend.position="bottom")
p3 <- p3a + p3b + plot_layout(widths=c(1,1.6))
ggsave(file.path(OUT,"Figure3_PGLS.pdf"), p3, width=11, height=5)
ggsave(file.path(OUT,"Figure3_PGLS.png"), p3, width=11, height=5, dpi=200)

## ---------- Figure 4: H1 x H3 synthesis (two axes) ----------
# ECM presence rate + mean per family from pf (ECM vs SAP+WD)
fam_lv <- c("GH7","GH6","GH28","GH5")
syn <- data.frame(family=fam_lv,
  K=c(0.46,0.80,0.87,1.00),
  ecm_pres=c(13,11,74,100),          # % ECM species with >=1 intact copy
  ecm_mean=c(0.22,0.21,2.12,11.97),  # mean intact copies in ECM
  h3p=c(8.7e-4,1.3e-4,2.9e-5,1.3e-7))
syn$family <- factor(syn$family, levels=fam_lv)
syn$type <- c("dedicated (lost + relaxed)","dedicated (lost + relaxed)","dedicated (intermediate)","multifunctional (kept + maintained)")
p4 <- ggplot(syn, aes(x=K, y=ecm_pres)) +
  annotate("rect", xmin=0.3, xmax=0.92, ymin=-3, ymax=40, alpha=0.06, fill="#D7263D") +
  annotate("rect", xmin=0.92, xmax=1.05, ymin=60, ymax=105, alpha=0.06, fill="#1B998B") +
  geom_vline(xintercept=1, linetype=2, color="grey60") +
  geom_point(aes(size=ecm_mean, color=K<1)) +
  geom_text_repel(aes(label=sprintf("%s (K=%.2f)", family, K)), size=3.4, box.padding=0.7, seed=1) +
  scale_color_manual(values=c(`TRUE`="#D7263D",`FALSE`="#1B998B"), guide="none") +
  scale_size_continuous(range=c(4,12), name="ECM mean\nintact copies") +
  scale_x_continuous(limits=c(0.35,1.05)) + scale_y_continuous(limits=c(-3,105)) +
  labs(x="H1: RELAX relaxation parameter K  (← more relaxed)",
       y="H3: % of ECM species retaining the family",
       title="Two axes of erosion: gene loss (H3) vs relaxed selection (H1)",
       subtitle="Dedicated cellulases GH6/GH7: lost AND relaxed  ·  multifunctional GH5: kept AND maintained") +
  th
ggsave(file.path(OUT,"Figure4_H1xH3_synthesis.pdf"), p4, width=8, height=6)
ggsave(file.path(OUT,"Figure4_H1xH3_synthesis.png"), p4, width=8, height=6, dpi=200)

cat("figures written to", OUT, "\n"); cat(list.files(OUT), sep="\n"); cat("\n")
