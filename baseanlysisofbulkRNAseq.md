# Bulk RNA-seq Analysis Pipeline (Part 1)

## 1. Bulk RNA-seq: Merge TPM Matrices and Annotate Protein-Coding Genes (Python)

```python
import os
import pandas as pd

# 1) Merge all RSEM .genes.results files in the current directory into a TPM matrix
files = [f for f in os.listdir('.') if f.endswith('.genes.results')]
print(f"Found {len(files)} files: {files}")

dfs = []
for f in files:
    sample_name = f.replace('.genes.results', '')
    df = pd.read_csv(f, sep='\t', usecols=['gene_id', 'TPM'])
    df.rename(columns={'TPM': sample_name}, inplace=True)
    dfs.append(df)

expr_matrix = dfs[0]
for df in dfs[1:]:
    expr_matrix = pd.merge(expr_matrix, df, on='gene_id', how='outer')

expr_matrix = expr_matrix.fillna(0)
expr_matrix.to_csv('All_samples_TPM_matrix.tsv', sep='\t', index=False)
print("✅ Merged TPM matrix saved as All_samples_TPM_matrix.tsv")

# 2) Keep protein-coding genes and collapse to gene symbols
annot_file = "/data_result/dengys/bulk/autobulk/ref/extracted_gene_info_mmu.csv"
output_file = "expr_matrix_geneSymbol_proteinCoding_sum.csv"

expr = expr_matrix
annot = pd.read_csv(annot_file)

expr.index = expr['gene_id'].astype(str)
annot["Gene_ID"] = annot["Gene_ID"].astype(str)

annot = annot[annot["Gene_Type"] == "protein_coding"]

expr = expr.merge(
    annot[["Gene_ID", "Gene_Name"]],
    left_index=True,
    right_on="Gene_ID",
    how="inner"
)

expr = expr.drop(columns=["Gene_ID"]).set_index("Gene_Name")
expr = expr.groupby(expr.index).sum()

expr.to_csv(output_file)
print(f"✅ Done! Result saved to: {output_file}")
print(f"Final matrix contains {expr.shape[0]} protein-coding genes")
```

---

## 2. Bulk RNA-seq: Sample-level t-SNE Visualization (R)

```r
## Dependencies
# install.packages(c("Rtsne","ggplot2","ggrepel","ggthemes"))
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(ggthemes)

## 0) Input matrix: rows = genes, columns = samples
# df <- read.csv("expected_count.genename_collapsed.filtered.csv", row.names = 1)
mat0 <- as.matrix(df)

## 1) CPM + log1p normalization + highly variable gene selection
libsize <- colSums(mat0)
cpm     <- t(t(mat0) / (libsize / 1e6))
logcpm  <- log1p(cpm)

gene_sd <- apply(logcpm, 1, sd, na.rm = TRUE)
mat_nzv <- logcpm[gene_sd > 0, , drop = FALSE]

N    <- min(2000, nrow(mat_nzv))
vars <- apply(mat_nzv, 1, var, na.rm = TRUE)
mat_hvg <- mat_nzv[order(vars, decreasing = TRUE)[1:N], , drop = FALSE]

## 2) PCA → t-SNE input
expr        <- t(mat_hvg)          # rows = samples
expr_scaled <- scale(expr)
pca_res     <- prcomp(expr_scaled, center = TRUE, scale. = FALSE)

K       <- min(50, ncol(pca_res$x))
tsne_in <- pca_res$x[, 1:K, drop = FALSE]

## 3) t-SNE
set.seed(123)
n     <- nrow(tsne_in)
perp  <- max(5, min(30, floor((n - 1) / 3)))
tsne_res <- Rtsne(tsne_in, dims = 2,
                  perplexity = perp,
                  verbose = TRUE,
                  check_duplicates = FALSE)

tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$Sample   <- rownames(expr)
tsne_df$Genotype <- ifelse(grepl("^Trem2_KO", tsne_df$Sample), "KO", "WT")
tsne_df$Time     <- sub(".*_(P\\d+)$", "\\1", tsne_df$Sample)

## 4) Plot
time_colors <- c("P7" = "#84a494",
                 "P14" = "#d4c464",
                 "P21" = "#f4a494")
geno_shapes <- c("WT" = 16, "KO" = 17)

p_tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2,
                              color = Time, shape = Genotype)) +
  geom_point(size = 4) +
  scale_color_manual(values = time_colors) +
  scale_shape_manual(values = geno_shapes) +
  theme_few(base_size = 14) +
  labs(
    title = sprintf("t-SNE (perplexity=%d, input=%d PCs, HVG=%d)",
                    perp, K, ncol(expr)),
    color = "Time point", shape = "Genotype"
  )

print(p_tsne)
ggsave("tsne_samples.png", p_tsne, width = 6, height = 5.2, dpi = 300)
```

---

## 3. Trem2-KO vs WT: Stage-wise DESeq2 Differential Expression (P7 / P14 / P21)

```r
library(DESeq2)

## Input: df is raw counts (genes x samples)
counts <- round(as.matrix(df))

samples <- colnames(counts)
colData <- data.frame(
  Sample   = samples,
  Genotype = ifelse(grepl("^Trem2_KO", samples), "KO", "WT"),
  Time     = sub(".*_(P\\d+)$", "\\1", samples),
  row.names = samples
)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData  = colData,
                              design   = ~ Genotype)
dds <- DESeq(dds)

## P7: KO vs WT
dds_P7  <- dds[, dds$Time == "P7"]
dds_P7  <- DESeq(dds_P7)
res_P7  <- results(dds_P7, contrast = c("Genotype", "KO", "WT"))
save(res_P7, file = "res_P7.rdata")

## P14: KO vs WT
dds_P14 <- dds[, dds$Time == "P14"]
dds_P14 <- DESeq(dds_P14)
res_P14 <- results(dds_P14, contrast = c("Genotype", "KO", "WT"))
save(res_P14, file = "res_P14.rdata")

## P21: KO vs WT
dds_P21 <- dds[, dds$Time == "P21"]
dds_P21 <- DESeq(dds_P21)
res_P21 <- results(dds_P21, contrast = c("Genotype", "KO", "WT"))
save(res_P21, file = "res_P21.rdata")

## Optional: write CSV
write.csv(as.data.frame(res_P7),  "res_P7.csv")
write.csv(as.data.frame(res_P14), "res_P14.csv")
write.csv(as.data.frame(res_P21), "res_P21.csv")
```

---

## 4. WT Developmental DEG and Stage Marker Definition (P7 / P14 / P21)

```r
library(DESeq2)
library(dplyr)

count <- read.csv("counts.csv", row.names = 1)

meta <- data.frame(
  sample   = colnames(count),
  genotype = ifelse(grepl("^WT", colnames(count)), "WT", "KO"),
  time     = sub(".*_P(\\d+)$", "P\\1", colnames(count))
)
rownames(meta) <- meta$sample

count_filtered <- round(count[rowSums(count) > 0, ])

## WT samples only, developmental design
dds <- DESeqDataSetFromMatrix(
  countData = count_filtered[, meta$genotype == "WT"],
  colData   = meta[meta$genotype == "WT", ],
  design    = ~ time
)
dds <- DESeq(dds)

res_P7   <- results(dds, contrast = c("time", "P7",  "P14"))
res_P7b  <- results(dds, contrast = c("time", "P7",  "P21"))
res_P14a <- results(dds, contrast = c("time", "P14", "P7"))
res_P14b <- results(dds, contrast = c("time", "P14", "P21"))
res_P21a <- results(dds, contrast = c("time", "P21", "P7"))
res_P21b <- results(dds, contrast = c("time", "P21", "P14"))

## Define stage markers via bidirectional contrasts
## (log2FC > 1, padj < 0.05, baseMean >= 5)
P7_markers <- intersect(
  rownames(res_P7[ which(res_P7$log2FoldChange  >  0.9 &
                          res_P7$pvalue         < 0.05 &
                          res_P7$baseMean     >= 5 ), ]),
  rownames(res_P7b[which(res_P7b$log2FoldChange > 0.9 &
                          res_P7b$pvalue        < 0.05 ), ])
)

P14_markers <- intersect(
  rownames(res_P14a[which(res_P14a$log2FoldChange > 0.9 &
                           res_P14a$pvalue        < 0.05 ), ]),
  rownames(res_P14b[which(res_P14b$log2FoldChange > 1 &
                           res_P14b$pvalue        < 0.05 ), ])
)

P21_markers <- intersect(
  rownames(res_P21a[which(res_P21a$log2FoldChange > 0.9 &
                           res_P21a$pvalue        < 0.05), ]),
  rownames(res_P21b[which(res_P21b$log2FoldChange > 1 &
                           res_P21b$pvalue         < 0.05  ), ])
)
```

---

## 5. Developmental Stage Ternary Plot (Based on WT Stage Markers)

```r
library(ggtern)
library(dplyr)
library(ggrepel)

data <- read.csv("expr_matrix_geneSymbol_proteinCoding_sum.csv",
                 header = TRUE, row.names = 1)
data <- data[, -1]                        # remove Ensembl ID column
data <- data[rowSums(data) > 0, ]

## Fix P7 sample naming if mis-labeled
colnames(data) <- gsub("WT([123])_P7",        "TMP\\1_P7",  colnames(data))
colnames(data) <- gsub("Trem2_KO([123])_P7",  "WT\\1_P7",   colnames(data))
colnames(data) <- gsub("TMP([123])_P7",       "Trem2_KO\\1_P7", colnames(data))

expr <- data

## 1) Make stage markers mutually exclusive
p7_only  <- setdiff(P7_markers,  union(P14_markers, P21_markers))
p14_only <- setdiff(P14_markers, union(P7_markers,  P21_markers))
p21_only <- setdiff(P21_markers, union(P7_markers,  P14_markers))

marker_all <- unique(c(p7_only, p14_only, p21_only))
marker_all <- intersect(marker_all, rownames(expr))

## 2) Row-wise z-score
expr_z <- t(scale(t(expr[marker_all, , drop = FALSE])))
expr_z[is.na(expr_z)] <- 0

## 3) Mean z-score per stage
score_z <- data.frame(
  sample = colnames(expr_z),
  P7  = if (length(p7_only)  > 0)
    colMeans(expr_z[rownames(expr_z) %in% p7_only,  , drop = FALSE]) else 0,
  P14 = if (length(p14_only) > 0)
    colMeans(expr_z[rownames(expr_z) %in% p14_only, , drop = FALSE]) else 0,
  P21 = if (length(p21_only) > 0)
    colMeans(expr_z[rownames(expr_z) %in% p21_only, , drop = FALSE]) else 0
)

## 4) Softmax → composition
softmax <- function(M, temp = 0.7){
  E <- exp(M / temp)
  E / rowSums(E)
}
prob <- as.data.frame(softmax(as.matrix(score_z[, c("P7","P14","P21")]), temp = 0.7))
score_tern <- cbind(score_z["sample"], prob)

score_tern$Genotype <- ifelse(grepl("^WT", score_tern$sample), "WT", "Trem2_KO")
score_tern$Time     <- gsub(".*_(P[0-9]+)$", "\\1", score_tern$sample)

score_tern$Genotype <- factor(score_tern$Genotype, levels = c("WT", "Trem2_KO"))
score_tern$Time     <- factor(score_tern$Time,     levels = c("P7", "P14", "P21"))

time_colors <- c("P7" = "#84a494", "P14" = "#d4c464", "P21" = "#f4a494")
geno_shapes <- c("WT" = 16, "Trem2_KO" = 17)

p <- ggtern(
  data = score_tern,
  aes(x = P14, y = P7, z = P21, color = Time, shape = Genotype)
) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = sample), size = 3, color = "black") +
  scale_color_manual(values = time_colors) +
  scale_shape_manual(values = geno_shapes) +
  theme_bw(base_size = 14) +
  theme_showarrows() +
  theme_nomask() +
  theme(
    legend.position       = "bottom",
    tern.axis.title.T     = element_text(size = 13, face = "bold"),
    tern.axis.title.L     = element_text(size = 13, face = "bold"),
    tern.axis.title.R     = element_text(size = 13, face = "bold"),
    plot.title            = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "Developmental Stage Composition",
    T = "P7",
    L = "P14",
    R = "P21",
    color = "Stage",
    shape = "Genotype"
  )

print(p)
```

---

## 6. KO vs WT (Three Stages) DEG Intersections and GO BP Enrichment (clusterProfiler)

```r
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

## Inputs: res_p7, res_p14, res_p21 (KO vs WT)
p_cut        <- 0.05
lfc_abs      <- 0.9
base_mean_min <- 5

filt_set <- function(res_df, dir = c("up","down"),
                     p_cut = 0.05,
                     lfc_abs = 0.9,
                     base_mean_min = 5){
  dir <- match.arg(dir)
  res_df <- as.data.frame(res_df)
  res_df$gene <- rownames(res_df)
  res_df <- res_df %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange),
           baseMean >= base_mean_min, pvalue < p_cut)
  if (dir == "up")   res_df <- res_df %>% filter(log2FoldChange >=  lfc_abs)
  if (dir == "down") res_df <- res_df %>% filter(log2FoldChange <= -lfc_abs)
  unique(res_df$gene)
}

bg_stage <- function(res_df, base_mean_min = 5){
  res_df <- as.data.frame(res_df)
  res_df$gene <- rownames(res_df)
  res_df %>%
    filter(!is.na(pvalue), baseMean >= base_mean_min) %>%
    pull(gene) %>%
    unique()
}

map2entrez <- function(genes){
  if (length(genes) == 0) return(character(0))
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Mm.eg.db) %>%
    pull(ENTREZID) %>%
    unique()
}

run_bp <- function(genes_entrez, universe_entrez){
  if (length(genes_entrez) == 0) return(NULL)
  enrichGO(gene          = genes_entrez,
           universe      = universe_entrez,
           OrgDb         = org.Mm.eg.db,
           ont           = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
}

safe_write <- function(obj, file){
  if (is.null(obj)) {
    write.csv(data.frame(), file, row.names = FALSE)
  } else {
    df <- as.data.frame(obj)
    write.csv(df, file, row.names = FALSE)
  }
}

## 1) Up- and down-regulated genes across three stages
up_P7    <- filt_set(res_p7,  "up",   p_cut, lfc_abs, base_mean_min)
up_P14   <- filt_set(res_p14, "up",   p_cut, lfc_abs, base_mean_min)
up_P21   <- filt_set(res_p21, "up",   p_cut, lfc_abs, base_mean_min)

down_P7  <- filt_set(res_p7,  "down", p_cut, lfc_abs, base_mean_min)
down_P14 <- filt_set(res_p14, "down", p_cut, lfc_abs, base_mean_min)
down_P21 <- filt_set(res_p21, "down", p_cut, lfc_abs, base_mean_min)

core_up   <- Reduce(intersect, list(up_P7, up_P14, up_P21))
core_down <- Reduce(intersect, list(down_P7, down_P14, down_P21))
ext_up    <- unique(c(up_P7, up_P14, up_P21))
ext_down  <- unique(c(down_P7, down_P14, down_P21))

cat("Core up:", length(core_up),
    " Core down:", length(core_down), "\n")
cat("Union up:", length(ext_up),
    " Union down:", length(ext_down), "\n")

## 2) Background = intersection of detectable genes
bg_p7  <- bg_stage(res_p7,  base_mean_min)
bg_p14 <- bg_stage(res_p14, base_mean_min)
bg_p21 <- bg_stage(res_p21, base_mean_min)
bg_all <- Reduce(intersect, list(bg_p7, bg_p14, bg_p21))

bg_entrez   <- map2entrez(bg_all)
core_up_e   <- map2entrez(core_up)
core_down_e <- map2entrez(core_down)
ext_up_e    <- map2entrez(ext_up)
ext_down_e  <- map2entrez(ext_down)

dir.create("GO_BP_ALL", showWarnings = FALSE)
writeLines(core_up,   "GO_BP_ALL/core_up_genes.txt")
writeLines(core_down, "GO_BP_ALL/core_down_genes.txt")
writeLines(ext_up,    "GO_BP_ALL/union_up_genes.txt")
writeLines(ext_down,  "GO_BP_ALL/union_down_genes.txt")

ego_core_up_bp   <- run_bp(core_up_e,   bg_entrez)
ego_core_down_bp <- run_bp(core_down_e, bg_entrez)
ego_ext_up_bp    <- run_bp(ext_up_e,    bg_entrez)
ego_ext_down_bp  <- run_bp(ext_down_e,  bg_entrez)

safe_write(ego_core_up_bp,   "GO_BP_ALL/GO_BP_core_up.csv")
safe_write(ego_core_down_bp, "GO_BP_ALL/GO_BP_core_down.csv")
safe_write(ego_ext_up_bp,    "GO_BP_ALL/GO_BP_union_up.csv")
safe_write(ego_ext_down_bp,  "GO_BP_ALL/GO_BP_union_down.csv")

cat("Four BP enrichment results written to directory: GO_BP_ALL/\n")
```

---

## 7. Global GO BP Dotplot (KO vs WT)

```r
library(ggplot2)
library(dplyr)
library(ggthemes)
library(stringr)
library(scales)

## go_bp is an enrichGO result (e.g., ego_core_up_bp, ego_ext_up_bp, etc.)
go_bp <- as.data.frame(GO)  # replace GO with the specific enrichGO object

go_bp$Description <- stringr::str_to_sentence(go_bp$Description)

to_ratio <- function(x){
  if (is.numeric(x)) return(x)
  sapply(strsplit(as.character(x), "/"),
         function(v) as.numeric(v[1]) / as.numeric(v[2]))
}
go_bp$GeneRatio_num <- to_ratio(go_bp$GeneRatio)
go_bp$negLog10Padj  <- -log10(go_bp$p.adjust + 1e-300)

TopN   <- 20
go_bp2 <- go_bp %>%
  arrange(desc(GeneRatio_num), p.adjust) %>%
  slice_head(n = TopN)

go_bp2$Description <- factor(go_bp2$Description,
                             levels = rev(go_bp2$Description))

palette6 <- c("#547298", "#8E9DBA", "#DADEE7",
              "#F1DDD6", "#DA9087", "#D24846")
val6     <- seq(0, 1, length.out = length(palette6))

pdf("GO_BP_ALL_KOvsWT.pdf", width = 10, height = 5)
ggplot(go_bp2, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(color = negLog10Padj, size = Count), alpha = 0.8) +
  scale_color_gradientn(
    colours = palette6,
    values  = val6,
    name    = "-log10(FDR)"
  ) +
  scale_size_continuous(range = c(3, 10), name = "Gene count") +
  labs(
    x = "Gene ratio",
    y = "GO term",
    title = "GO enrichment (BP)"
  ) +
  theme_few() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )
dev.off()
```

---

## 8. Number of DEGs per Developmental Stage (P7 / P14 / P21)

```r
library(ggplot2)
library(ggthemes)

df <- data.frame(
  Stage = rep(c("P7", "P14", "P21"), each = 2),
  Regulation = rep(c("Up", "Down"), times = 3),
  Count = c(71, 858,   # P7
            699, 231,  # P14
            26, 743)   # P21
)

df$Stage      <- factor(df$Stage,      levels = c("P7", "P14", "P21"))
df$Regulation <- factor(df$Regulation, levels = c("Up", "Down"))

colors <- c("Up" = "#DA9087", "Down" = "#547298")

p <- ggplot(df, aes(x = Stage, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Developmental Stage",
       y = "Number of DEGs",
       fill = "Regulation",
       title = "Up- and Down-regulated Genes") +
  theme_few() +
  theme(
    axis.title.x  = element_text(size = 16, vjust = -0.5),
    axis.title.y  = element_text(size = 16, vjust = 1.5),
    axis.text.x   = element_text(size = 13),
    axis.text.y   = element_text(size = 12),
    plot.title    = element_text(size = 17, hjust = 0.5),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, size = 1)
  )

ggsave("DEG_barplot.png", p, width = 7, height = 6, dpi = 300)
```

---

## 9. GSVA: GO BP Pathway Scores and Heatmaps

```r
library(GSEABase)
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

Go_bp_path <- "D:/GSVA/msigdb_v2025.1.Mm_files_to_download_locally/" %+%
              "msigdb_v2025.1.Mm_GMTs/m5.go.bp.v2025.1.Mm.symbols.gmt"
geneset <- getGmt(Go_bp_path)

data <- read.csv("expr_matrix_geneSymbol_proteinCoding_sum.csv",
                 header = TRUE, row.names = 1)
data <- data[, -1]                       # remove Ensembl column

## Fix P7 naming
colnames(data) <- gsub("WT([123])_P7",        "TMP\\1_P7",  colnames(data))
colnames(data) <- gsub("Trem2_KO([123])_P7",  "WT\\1_P7",   colnames(data))
colnames(data) <- gsub("TMP([123])_P7",       "Trem2_KO\\1_P7", colnames(data))

data <- data[rowSums(data) > 0, ]
data <- na.omit(data)

expr_mat <- as.matrix(data)
expr_mat <- log2(expr_mat + 1)

## Filter gene sets
feat <- rownames(expr_mat)
gs_list <- geneIds(geneset)
gs_list_filtered <- lapply(gs_list, function(v) intersect(v, feat))
gs_list_filtered <- gs_list_filtered[
  sapply(gs_list_filtered, function(v) length(v) >= 5 & length(v) <= 5000)
]
cat("Number of retained gene sets:", length(gs_list_filtered), "\n")

param <- gsvaParam(
  exprData = expr_mat,
  geneSets = gs_list_filtered,
  kcdf     = "Gaussian",
  minSize  = 5,
  maxSize  = 5000
)
gsva_result <- gsva(param)

## Column annotation
samples  <- colnames(gsva_result)
stage    <- ifelse(grepl("P7",  samples), "P7",
                   ifelse(grepl("P14", samples), "P14", "P21"))
genotype <- ifelse(grepl("^WT", samples), "WT", "KO")

ann_col <- HeatmapAnnotation(
  Stage = stage,
  Genotype = genotype,
  col = list(
    Stage    = c(P21 = "#D39394", P14 = "#E5CC96", P7 = "#84a494"),
    Genotype = c(WT  = "#DA9087", KO  = "#547298")
  ),
  annotation_name_side = "left",
  gp = gpar(col = NA)
)

## Select variable pathways and row-wise z-score
topN <- 100
keep <- head(order(rowVars(gsva_result), decreasing = TRUE), topN)
mat  <- gsva_result[keep, , drop = FALSE]
mat_z <- t(scale(t(mat)))

split_fac <- factor(stage, levels = c("P7","P14","P21"))
ord <- order(
  factor(stage,    levels = c("P7", "P14", "P21")),
  factor(genotype, levels = c("WT", "KO"))
)

mat_z     <- mat_z[, ord]
split_fac <- split_fac[ord]
ann_col   <- ann_col[ord]

col_fun <- colorRamp2(
  c(-1.5, -0.9, -0.3, 0.3, 0.9, 1.5),
  c("#547298", "#8E9DBA", "#DADEE7",
    "#F1DDD6", "#DA9087", "#D24846")
)

## Heatmap with legend only
ht_with_legend <- Heatmap(
  mat_z,
  name = "GSVA z",
  col = col_fun,
  top_annotation   = ann_col,
  column_split     = split_fac,
  cluster_columns  = FALSE,
  cluster_rows     = TRUE,
  show_column_names = FALSE,
  show_row_names   = FALSE,
  border           = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    direction = "vertical",
    title_position = "topcenter",
    legend_height   = unit(4, "cm")
  )
)

pdf("gsva_heatmap_with_legend.pdf", width = 10, height = 8)
draw(ht_with_legend,
     heatmap_legend_side     = "right",
     annotation_legend_side  = "right")
dev.off()

## Heatmap with sparse row labels
ht_temp <- Heatmap(
  mat_z,
  name = "GSVA z",
  col = col_fun,
  top_annotation  = ann_col,
  column_split    = split_fac,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  show_row_names  = FALSE,
  border          = FALSE
)
ht_temp   <- draw(ht_temp)
row_order <- row_order(ht_temp)

show_positions <- seq(1, length(row_order), by = 3)
row_labels     <- rep("", nrow(mat_z))
row_labels[row_order[show_positions]] <- rownames(mat_z)[row_order[show_positions]]

ht_with_labels <- Heatmap(
  mat_z,
  name = "GSVA z",
  col = col_fun,
  top_annotation   = ann_col,
  column_split     = split_fac,
  cluster_columns  = FALSE,
  cluster_rows     = TRUE,
  show_column_names = FALSE,
  show_row_names   = TRUE,
  row_labels       = row_labels,
  row_names_gp     = gpar(fontsize = 7),
  row_names_max_width = unit(12, "cm"),
  border           = FALSE,
  show_heatmap_legend = FALSE
)

pdf("gsva_heatmap_with_labels.pdf", width = 12, height = 12)
draw(ht_with_labels,
     show_annotation_legend = FALSE,
     padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()
```

---

## 10. P14 / P21 Volcano Plot and GO/KEGG Enrichment (Example: P14)

```r
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggthemes)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Hmisc)
library(scales)

## 1) P14: Volcano plot
res_df <- as.data.frame(res_p14)
res_df$gene <- rownames(res_df)

min_pval <- 1e-300

res_df <- res_df %>%
  filter(!is.na(pvalue),
         !is.na(log2FoldChange),
         !is.na(baseMean)) %>%
  mutate(
    pvalue_adj = pmax(pvalue, min_pval),
    log10p     = -log10(pvalue_adj),
    direction  = case_when(
      log2FoldChange >=  0.9 & pvalue <= 0.05 ~ "Up",
      log2FoldChange <= -0.9 & pvalue <= 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

up_count   <- sum(res_df$direction == "Up",   na.rm = TRUE)
down_count <- sum(res_df$direction == "Down", na.rm = TRUE)

label_genes_up <- res_df %>%
  filter(pvalue < 0.05,
         log2FoldChange > 2,
         baseMean > 300,
         !grepl("^ens|^gm|^rps|^rpl", gene, ignore.case = TRUE),
         !grepl("rik$", gene, ignore.case = TRUE)) %>%
  mutate(rank_score = 4 * rank(-abs(log2FoldChange)) +
                      0.5 * rank(-log10(pvalue_adj))) %>%
  arrange(rank_score) %>%
  slice_head(n = 10)

label_genes_down <- res_df %>%
  filter(pvalue < 0.05,
         log2FoldChange < -1,
         baseMean > 500,
         !grepl("^ens|^gm|^rps|^rpl", gene, ignore.case = TRUE),
         !grepl("rik$", gene, ignore.case = TRUE)) %>%
  mutate(rank_score = 4 * rank(-abs(log2FoldChange)) +
                      0.5 * rank(-log10(pvalue_adj))) %>%
  arrange(rank_score) %>%
  slice_head(n = 10)

label_genes <- bind_rows(label_genes_up, label_genes_down)

p <- ggplot(res_df, aes(x = log2FoldChange, y = log10p)) +
  geom_point(aes(color = log2FoldChange, size = baseMean), alpha = 0.7) +
  geom_vline(xintercept = c(-0.9, 0.9), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_gradientn(
    colours = colorRampPalette(c(
      "#547298", "#8E9DBA", "#DADEE7",
      "#F1DDD6", "#DA9087", "#D24846"
    ))(100),
    values = rescale(c(-5, -3, -1, 0, 1, 3, 5)),
    limits = c(-5, 5),
    oob    = scales::squish,
    name   = expression(log[2]*"FC")
  ) +
  scale_size_continuous(
    range = c(0.5, 6),
    name  = "Base Mean",
    trans = "log10",
    breaks = c(10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 4.5,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0
  ) +
  annotate("text", x = -5, y = 250,
           label = paste0("Down (", down_count, ")"),
           hjust = 0, color = "#547298", size = 6) +
  annotate("text", x = 3.5, y = 250,
           label = paste0("Up (", up_count, ")"),
           hjust = 0, color = "#DA9087", size = 6) +
  labs(
    x = expression(log[2]*"FoldChange (KO vs WT)"),
    y = expression(-log[10]*italic(P)*"-value"),
    title = "P14: KO vs WT"
  ) +
  theme_few() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "right",
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black"),
    axis.text     = element_text(size = 16, color = "black"),
    axis.title    = element_text(face = "bold", size = 18)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 310)) +
  scale_x_continuous(limits = c(-5, 5))

pdf("Volcano_P14_KOvsWT.pdf", width = 8, height = 6)
print(p)
dev.off()

## 2) P14 upregulated genes GO BP enrichment
up <- subset(res_p14,
             subset = log2FoldChange >= 0.9 & pvalue <= 0.05)
genes_up        <- rownames(up)
genes_up_entrez <- bitr(genes_up,
                        fromType = "SYMBOL",
                        toType   = "ENTREZID",
                        OrgDb    = org.Mm.eg.db)

go_bp <- enrichGO(
  gene          = genes_up_entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_bp <- as.data.frame(go_bp)
go_bp$Description <- capitalize(go_bp$Description)

go_bp2 <- head(go_bp, 10) %>%
  arrange(desc(GeneRatio))

go_bp2$GeneRatio <- as.numeric(
  sapply(strsplit(as.character(go_bp2$GeneRatio), "/"),
         function(x) as.numeric(x[1]) / as.numeric(x[2]))
)
go_bp2$Description <- factor(go_bp2$Description,
                             levels = rev(go_bp2$Description))

max_padj <- max(go_bp2$p.adjust, na.rm = TRUE)
min_padj <- min(go_bp2$p.adjust, na.rm = TRUE)

pdf("GO_BP_P14_KOvsWT.pdf", width = 8, height = 6)
ggplot(go_bp2, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(color = p.adjust, size = Count), alpha = 0.7) +
  scale_color_gradientn(
    colours = rev(c(
      "#547298", "#8E9DBA", "#DADEE7",
      "#F1DDD6", "#DA9087", "#D24846"
    )),
    values = rescale(c(min_padj,
                       min_padj * 5,
                       max_padj * 0.7,
                       max_padj * 1.2)),
    limits = c(min_padj, max_padj),
    oob    = scales::squish,
    name   = "Adjusted P-value"
  ) +
  scale_size_continuous(
    range = c(3, 10), name = "Gene Count"
  ) +
  theme_few() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  ) +
  labs(
    x = "Gene Ratio",
    y = "GO Term",
    title = "GO Enrichment Analysis"
  )
dev.off()

## 3) P14 upregulated genes KEGG enrichment and barplot
kegg_enrich <- enrichKEGG(
  gene         = genes_up_entrez$ENTREZID,
  organism     = "mmu",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

go_data <- as.data.frame(kegg_enrich)
go_data$log10_adjustp <- -log10(go_data$p.adjust)
go_data$Description   <- capitalize(go_data$Description)

go_data2 <- head(go_data[order(go_data$p.adjust), ], 10)
go_data2$GeneRatio <- as.numeric(
  sapply(strsplit(as.character(go_data2$GeneRatio), "/"),
         function(x) as.numeric(x[1]) / as.numeric(x[2]))
)
go_data2 <- go_data2 %>% arrange(GeneRatio)
go_data2$Description <- factor(go_data2$Description,
                               levels = go_data2$Description)

pdf("KEGG_enrichment_plot.pdf", width = 7, height = 5)
ggplot(go_data2, aes(x = GeneRatio, y = Description)) +
  geom_bar(stat = "identity",
           aes(fill = log10_adjustp),
           width = 0.7) +
  scale_fill_gradientn(
    colours = c(
      "#547298", "#8E9DBA", "#DADEE7",
      "#F1DDD6", "#DA9087", "#D24846"
    ),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
    limits = c(0, 19),
    oob    = scales::squish,
    name   = expression(-log[10]*"(Adjusted P-value)")
  ) +
  geom_text(aes(x = 0.01, label = Description),
            hjust = 0, color = "black", size = 5) +
  labs(x = "Gene Ratio", y = "") +
  theme_minimal() +
  ggtitle("Top 10 KEGG Pathways") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    axis.line.x  = element_line(color = "black")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
dev.off()
```
