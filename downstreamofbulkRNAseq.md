# Bulk RNA-seq Analysis Pipeline (Part 2)

## 11. Bulk RNA-seq: Preprocessing for WGCNA

```r
## Load dependencies
library(WGCNA)
library(tidyverse)

allowWGCNAThreads()

## Input: normalized expression matrix (e.g. variance-stabilized or log2 CPM)
## rows = genes, columns = samples
expr <- read.csv("expr_matrix_geneSymbol_proteinCoding_sum.csv",
                 header = TRUE, row.names = 1, check.names = FALSE)

## Drop Ensembl ID column if present
if ("Ensembl" %in% colnames(expr)) {
  expr <- expr[, !(colnames(expr) %in% "Ensembl"), drop = FALSE]
}

## Keep genes with sufficient variability
expr <- as.data.frame(expr)
expr <- expr[rowSums(expr > 0) >= 3, ]
var_cutoff <- quantile(apply(expr, 1, var), 0.25)
expr_filt  <- expr[apply(expr, 1, var) > var_cutoff, ]

## Transpose for WGCNA (samples x genes)
datExpr <- t(expr_filt)
gsg     <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

## Sample dendrogram to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("WGCNA_sampleClustering.pdf", width = 8, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers",
     xlab = "", sub = "", cex = 0.7)
dev.off()
```

---

## 12. WGCNA: Network Construction and Module Detection

```r
## Build sample trait data frame
samples <- rownames(datExpr)

trait <- data.frame(
  Sample = samples,
  Genotype = ifelse(grepl("^WT",  samples), "WT", "KO"),
  Time     = gsub(".*_(P[0-9]+)$", "\\1", samples),
  stringsAsFactors = FALSE
)
trait$Group <- paste(trait$Genotype, trait$Time, sep = "_")
rownames(trait) <- trait$Sample

## Choose soft-thresholding power
powers  <- c(c(1:10), seq(12, 30, 2))
sft     <- pickSoftThreshold(datExpr, powerVector = powers,
                             networkType = "signed", verbose = 5)

pdf("WGCNA_pickSoftThreshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.8, col = "red")
abline(h = 0.8, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.8, col = "red")
dev.off()

## Choose a power that reaches R^2 ~ 0.8, or use a fixed value (e.g. 12)
softPower <- 12

adjacency <- adjacency(datExpr, power = softPower, type = "signed")
TOM      <- TOMsimilarity(adjacency)
dissTOM  <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf("WGCNA_geneClustering.pdf", width = 8, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

## Dynamic tree cut to define modules
minModuleSize  <- 30
dynamicMods    <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
dynamicColors  <- labels2colors(dynamicMods)

pdf("WGCNA_geneTree_dynamicModules.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## Calculate module eigengenes
MEs0 <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes
MEs  <- orderMEs(MEs0)
moduleColors <- dynamicColors
```

---

## 13. WGCNA: Module–Trait Relationships and Heatmap

```r
## One-hot encoding of 6 groups: WT/KO x P7/P14/P21
trait$Group <- factor(trait$Group,
                      levels = c("WT_P7", "KO_P7",
                                 "WT_P14","KO_P14",
                                 "WT_P21","KO_P21"))
designMat <- model.matrix(~ 0 + Group, data = trait)
colnames(designMat) <- levels(trait$Group)

## Correlation between module eigengenes and traits
moduleTraitCor    <- cor(MEs, designMat, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

## Build text matrix: correlation (p-value)
textMat <- paste0(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 2), ")")
dim(textMat) <- dim(moduleTraitCor)

pdf("WGCNA_moduleTraitRelationships.pdf", width = 10, height = 8)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(designMat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMat,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Module–trait relationships")
dev.off()
```

---

## 14. WGCNA: Extract Module Genes and Perform GO / KEGG Enrichment

```r
library(clusterProfiler)
library(org.Mm.eg.db)

## Example: select a biologically relevant module (e.g. 'turquoise')
moduleOfInterest <- "turquoise"

moduleGenes <- names(datExpr)[moduleColors == moduleOfInterest]

## Map to Entrez IDs
gene_df <- bitr(moduleGenes,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Mm.eg.db)

## Define background as all WGCNA genes
bg_df <- bitr(names(datExpr),
              fromType = "SYMBOL",
              toType   = "ENTREZID",
              OrgDb    = org.Mm.eg.db)

## GO BP enrichment
ego_bp <- enrichGO(gene          = unique(gene_df$ENTREZID),
                   universe      = unique(bg_df$ENTREZID),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

go_bp <- as.data.frame(ego_bp)
write.csv(go_bp, file = paste0("WGCNA_", moduleOfInterest, "_GO_BP.csv"),
          row.names = FALSE)

## KEGG enrichment
kegg_res <- enrichKEGG(gene         = unique(gene_df$ENTREZID),
                       universe     = unique(bg_df$ENTREZID),
                       organism     = "mmu",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

kegg_df <- as.data.frame(kegg_res)
write.csv(kegg_df,
          file = paste0("WGCNA_", moduleOfInterest, "_KEGG.csv"),
          row.names = FALSE)

## Simple KEGG barplot
if (nrow(kegg_df) > 0) {
  kegg_top <- head(kegg_df[order(kegg_df$p.adjust), ], 10)
  kegg_top$GeneRatio <- sapply(strsplit(as.character(kegg_top$GeneRatio), "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]))
  kegg_top$Description <- factor(kegg_top$Description,
                                 levels = rev(kegg_top$Description))

  pdf(paste0("WGCNA_", moduleOfInterest, "_KEGG_barplot.pdf"),
      width = 7, height = 5)
  ggplot(kegg_top, aes(x = GeneRatio, y = Description)) +
    geom_bar(stat = "identity", aes(fill = -log10(p.adjust)), width = 0.7) +
    scale_fill_gradientn(
      colours = c("#547298", "#8E9DBA", "#DADEE7",
                  "#F1DDD6", "#DA9087", "#D24846"),
      name = expression(-log[10]*"(Adjusted P-value)")
    ) +
    labs(x = "Gene ratio", y = "", title = "Top KEGG pathways") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x  = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10)
    )
  dev.off()
}
```

---

## 15. Time-series Clustering with Mfuzz (Group-averaged Expression)

```r
library(Mfuzz)
library(Biobase)
library(tidyverse)

## Expression matrix: genes x samples
expr <- read.csv("expr_matrix_geneSymbol_proteinCoding_sum.csv",
                 header = TRUE, row.names = 1, check.names = FALSE)
expr <- expr[, !(colnames(expr) %in% "Ensembl"), drop = FALSE]
expr <- expr[rowSums(expr) > 0, ]
expr_mat <- as.matrix(expr)

## Build group labels: WT/KO x P7/P14/P21
samples <- colnames(expr_mat)
Genotype <- ifelse(grepl("^WT", samples), "WT", "KO")
Time     <- gsub(".*_(P[0-9]+)$", "\\1", samples)
Group    <- paste(Genotype, Time, sep = "_")

## Group-averaged expression (genes x 6 groups)
groups <- unique(Group)
groups <- factor(groups, levels = c("WT_P7", "KO_P7",
                                    "WT_P14","KO_P14",
                                    "WT_P21","KO_P21"))

group_mat <- sapply(levels(groups), function(g) {
  colMeans(expr_mat[, Group == g, drop = FALSE], na.rm = TRUE)
})
colnames(group_mat) <- levels(groups)

## Build ExpressionSet for Mfuzz
eset <- ExpressionSet(assayData = as.matrix(group_mat))

## Standardize
eset_s <- standardise(eset)

## Estimate fuzzification parameter m
m <- mestimate(eset_s)

## Choose number of clusters, e.g. 6
set.seed(123)
k <- 6
cl <- mfuzz(eset_s, c = k, m = m)

## Plot cluster trajectories
pdf("Mfuzz_clusters.pdf", width = 10, height = 8)
mfuzz.plot(eset_s, cl = cl,
           mfrow = c(2, 3),
           time.labels = colnames(group_mat),
           new.window = FALSE,
           xlab = "Group (WT/KO x Stage)",
           ylab = "Standardised expression")
dev.off()

## Export cluster membership
cluster_membership <- as.data.frame(cl$cluster)
colnames(cluster_membership) <- "Cluster"
cluster_membership$Gene <- rownames(cluster_membership)
write.csv(cluster_membership,
          "Mfuzz_cluster_membership.csv", row.names = FALSE)
```

---

## 16. WT vs KO Expression Trends per Mfuzz Cluster

```r
library(tidyverse)

## Load Mfuzz membership and group-averaged matrix
cluster_membership <- read.csv("Mfuzz_cluster_membership.csv")
group_mat <- group_mat  # from previous section (genes x 6 groups)

## Long-format mean expression per cluster x group
expr_long <- as.data.frame(group_mat)
expr_long$Gene   <- rownames(expr_long)
expr_long        <- expr_long %>% left_join(cluster_membership, by = "Gene")

expr_long <- expr_long %>%
  pivot_longer(cols = starts_with("WT_") | starts_with("KO_"),
               names_to = "Group", values_to = "Expression") %>%
  separate(Group, into = c("Genotype","Time"), sep = "_", remove = FALSE) %>%
  mutate(Time_num = as.numeric(gsub("P", "", Time)))

## Summarize mean ± SE per cluster, genotype, time
expr_summary <- expr_long %>%
  group_by(Cluster, Genotype, Time_num) %>%
  summarise(
    Mean = mean(Expression, na.rm = TRUE),
    SE   = sd(Expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

## Plot
expr_summary$Time_num <- factor(expr_summary$Time_num,
                                levels = c(7, 14, 21),
                                labels = c("P7", "P14", "P21"))

pdf("WT_KO_Mfuzz_trends.pdf", width = 12, height = 8)
ggplot(expr_summary,
       aes(x = Time_num, y = Mean,
           color = Genotype, group = Genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE,
                  fill = Genotype, group = Genotype),
              alpha = 0.15, color = NA) +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_color_manual(values = c("WT" = "#DA9087", "KO" = "#547298")) +
  scale_fill_manual(values = c("WT" = "#DA9087", "KO" = "#547298")) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Developmental stage",
       y = "Mean expression",
       title = "WT vs KO expression trends by Mfuzz cluster")
dev.off()

write.csv(expr_summary,
          "WT_KO_cluster_mean_expression.csv",
          row.names = FALSE)
```

---

## 17. UMAP of Samples Colored by Group or Mfuzz-derived Signatures

```r
library(uwot)
library(ggplot2)
library(ggthemes)
library(tidyverse)

## Use the same group-averaged or sample-level matrix
expr_mat <- as.matrix(expr_filt)  # or expr_mat from previous sections

## Z-score per gene
expr_z <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0

## Run UMAP on samples
set.seed(123)
umap_res <- umap(t(expr_z), n_neighbors = 6, min_dist = 0.3, metric = "euclidean")
umap_df  <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample   <- colnames(expr_mat)
umap_df$Genotype <- ifelse(grepl("^WT", umap_df$Sample), "WT", "KO")
umap_df$Time     <- gsub(".*_(P[0-9]+)$", "\\1", umap_df$Sample)
umap_df$Group    <- paste(umap_df$Genotype, umap_df$Time, sep = "_")

## UMAP colored by Group
pal_group <- c("WT_P7"  = "#D4E4D9",
               "KO_P7"  = "#84A494",
               "WT_P14" = "#F0E1B9",
               "KO_P14" = "#D4C464",
               "WT_P21" = "#F6C4B4",
               "KO_P21" = "#F4A494")

umap_df$Group <- factor(umap_df$Group,
                        levels = c("WT_P7","KO_P7",
                                   "WT_P14","KO_P14",
                                   "WT_P21","KO_P21"))

p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Group)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = pal_group, drop = FALSE) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    legend.key      = element_blank()
  ) +
  labs(title = "UMAP of bulk samples colored by group")

pdf("UMAP_samples_by_group.pdf", width = 7, height = 6)
print(p1)
dev.off()
```
