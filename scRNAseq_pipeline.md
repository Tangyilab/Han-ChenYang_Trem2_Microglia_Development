# Microglia Development scRNA-seq Analysis Pipeline

This document summarizes the **single-cell RNA-seq** analysis workflow for microglial development, including:

1. Seurat object construction, QC, and merging across samples  
2. Doublet removal and integrated analysis (LogNormalize + scVI integration)  
3. Seurat → Monocle3 conversion and pseudotime analysis  

The code below is organized as three main scripts:

- `Seurat_QC_Merge_Main.R` – per-sample QC and merged raw object construction  
- `Seurat_Integration_Main.R` – doublet removal and multi-batch integration  
- `Seurat_to_Monocle3.R` – pseudotime analysis based on the integrated Seurat object  

---

## 1. Seurat QC and Merging Across Samples (`Seurat_QC_Merge_Main.R`)

```r
#### Seurat_QC_Merge_Main.R
### 主流程：逐样本构建 Seurat 对象并统一 QC，合并为 total_raw / total_flited

library(Seurat)
library(DropletQC)  # 用于 nuclear_fraction_tags 等指标
library(dplyr)

source("/data_result/dengys/scRNA/AutoScRNA/Script/R/create_seurat_with_qc.R")

basepath <- "/data/dengys/microglia_develop"
setwd(basepath)
dir.create("./Results/rawdata", recursive = TRUE)

## 读取样本信息
sampleinfo <- read.csv("sample_info.csv", header = TRUE, row.names = 1)
sampleid_list <- rownames(sampleinfo)

## 循环构建 Seurat 对象
len_sample <- nrow(sampleinfo)
scRNAlist <- list()

for (i in 1:len_sample) {
  # create_seurat_with_qc 内部完成：从 CellRanger / CellBender 结果导入、基础过滤和 QC 指标计算
  scRNA <- create_seurat_with_qc(sampleid_list[i], basepath)
  scRNAlist[[i]] <- scRNA
  print(i)
  rm(scRNA)
}

#### 合并所有样本
total_raw <- merge(scRNAlist[[1]], y = scRNAlist[2:len_sample])
Idents(total_raw) <- "orig.ident"

#### 打上 meta 信息（与 sample_info 表格对齐）
sampleinfo$sample <- row.names(sampleinfo)
meta_seu <- total_raw@meta.data
meta_seu$barcode <- row.names(meta_seu)
meta_seu <- merge(meta_seu, sampleinfo, by = "sample")
row.names(meta_seu) <- meta_seu$barcode
total_raw <- AddMetaData(total_raw, metadata = meta_seu)

### QC 前后 violin plot
pdf("./Results/rawdata/sampleinfo_raw.pdf", width = 16, height = 12)
VlnPlot(
  total_raw,
  features = c(
    "nFeature_RNA", "nCount_RNA",
    "nuclear_fraction", "background_fraction",
    "cell_probability", "cell_size", "droplet_efficiency",
    "percent.mt", "percent.hb"
  ),
  pt.size = 0,
  ncol = 3
)
dev.off()

## 将 layers 合并（Seurat v5）
total_raw <- JoinLayers(total_raw)
save(total_raw, file = "./Results/rawdata/total_raw.Rdata")

#### 过滤细胞（严格 QC 阈值）
total_flited <- subset(
  total_raw,
  subset =
    nFeature_RNA > 300 &
    nCount_RNA   > 1000 &
    nCount_RNA   < 10000 &
    percent.mt   < 10 &
    cell_probability > 0.5 &
    percent.hb   < 0.1 &
    cell_status  == "cell"
)

pdf("./Results/rawdata/sampleinfo_afterflited.pdf", width = 16, height = 12)
VlnPlot(
  total_flited,
  features = c(
    "nFeature_RNA", "nCount_RNA",
    "nuclear_fraction", "background_fraction",
    "cell_probability", "cell_size", "droplet_efficiency",
    "percent.mt", "percent.hb"
  ),
  pt.size = 0,
  ncol = 3
)
dev.off()

save(total_flited, file = "./Results/rawdata/total_flited.Rdata")
```

---

## 2. Doublet Removal and Multi-batch Integration (`Seurat_Integration_Main.R`)

```r
#### Seurat_Integration_Main.R
## 目标：加载多批次 total_flited 对象，去除双细胞，按 orig.ident 分批，使用 LogNormalize + scVI 实现跨批次整合

library(Seurat)
library(patchwork)
library(DropletUtils)
library(scDblFinder)
library(SeuratWrappers)
library(SingleCellExperiment)

### 加载必要函数
source("/data_result/dengys/scRNA/AutoScRNA/Script/R/filter_doublets.R")
source("/data_result/dengys/scRNA/AutoScRNA/Script/R/seurat_integration_pipeline.R")

basepath <- "/data/dengys/microglia_develop"
setwd(basepath)
dir.create("./Results/integrated_data", recursive = TRUE)

# 加载两批次的 total_flited 对象
raw_dir <- "./Results/rawdata/total_flited.Rdata"

# 批次 2（例如 late 或额外批次）
load("/data/dengys/microglia_develop/batch2/fastq_download/Results/rawdata/total_flited.Rdata")
total_flited_li <- total_flited
total_flited_li <- JoinLayers(total_flited_li)

# 批次 1（主批次）
load(raw_dir)
total_flited_ha <- total_flited
total_flited_ha <- JoinLayers(total_flited_ha)

### 合并两批数据
total_flited <- merge(total_flited_ha, total_flited_li)
total_flited <- JoinLayers(total_flited)

### 去除双细胞（自定义函数封装 scDblFinder / Scrublet 等）
total_nodouble <- no_double(total_flited)

###### 数据整合
## 将 RNA assay 按 orig.ident 划分为 list，适配 Seurat integration 接口
total_nodouble[["RNA"]] <- split(total_nodouble[["RNA"]], f = total_nodouble$orig.ident)

### 自定义整合方法
# norm_methods <- c("LogNormalize", "SCT")
# integrate_methods <- c("CCA", "RPCA", "Harmony", "scVI")
norm_methods      <- c("LogNormalize")
integrate_methods <- c("scVI")

## 执行整合（封装 PCA / neighbors / UMAP / clustering 等标准流程）
process_seurat_integrate(
  total_nodouble,
  norm_methods      = norm_methods,
  integrate_methods = integrate_methods,
  output_dir        = "./Results/integrated_data/",
  verbose           = TRUE
)

# 输出对象中通常包含：
# - integrated assay / latent scVI embedding
# - UMAP reductions (e.g. umap.integrated.scvi)
# - 统一的 cluster / annotation 元数据
```

---

## 3. Seurat → Monocle3 Pseudotime Analysis (`Seurat_to_Monocle3.R`)

```r
## ============================================================
## 🧬 Seurat v5 → Monocle3 转换 + 拟时序分析框架
## ============================================================

library(Seurat)
library(monocle3)
library(Matrix)
library(dplyr)

# 假设你的整合后 Seurat 对象叫：
# clustered_obj
# reductions: pca, integrated.scvi, umap.integrated.scvi

# --------------------------
# 1️⃣ 提取表达矩阵
# --------------------------
# 建议使用 counts（更适合拟时序模型），也可以根据需要改为 data slot
expr_mat <- clustered_obj@assays$RNA@counts
# 如果希望使用 log-normalized data：
# expr_mat <- GetAssayData(clustered_obj, assay = "RNA", slot = "data")

# --------------------------
# 2️⃣ 提取元数据
# --------------------------
cell_metadata <- clustered_obj@meta.data

# --------------------------
# 3️⃣ 提取特征信息（基因注释）
# --------------------------
gene_annotation <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names       = rownames(expr_mat)
)

# --------------------------
# 4️⃣ 构建 Monocle3 的 CDS 对象
# --------------------------
cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata   = cell_metadata,
  gene_metadata   = gene_annotation
)

# --------------------------
# 5️⃣ 导入 Seurat 的 UMAP 坐标（保持视图一致）
# --------------------------
if ("umap.integrated.scvi" %in% names(clustered_obj@reductions)) {
  umap_coords <- clustered_obj@reductions$umap.integrated.scvi@cell.embeddings
} else if ("umap" %in% names(clustered_obj@reductions)) {
  umap_coords <- clustered_obj@reductions$umap@cell.embeddings
} else {
  stop("UMAP reduction not found in Seurat object.")
}

# 将 UMAP 坐标传递到 cds 对象中
reducedDims(cds)$UMAP <- umap_coords

# --------------------------
# 6️⃣ 将 Seurat 分群标签和其它 meta 信息写入 CDS
# --------------------------
if ("seurat_clusters" %in% colnames(clustered_obj@meta.data)) {
  colData(cds)$seurat_clusters <- clustered_obj@meta.data$seurat_clusters
}

# 例如还可以写入 Genotype / Stage 等
# colData(cds)$Genotype <- clustered_obj@meta.data$Genotype
# colData(cds)$Stage    <- clustered_obj@meta.data$Stage

# --------------------------
# 7️⃣ Monocle3 图结构学习 + 伪时间计算
# --------------------------
cds <- learn_graph(cds)
cds <- order_cells(cds)  # 交互式 / 手动指定 root cells 时可在这里调整

# --------------------------
# 8️⃣ 可视化：UMAP 轨迹和伪时间
# --------------------------
# 伪时间着色
plot_cells(
  cds,
  color_cells_by         = "pseudotime",
  show_trajectory_graph  = TRUE,
  label_groups_by_cluster = FALSE,
  label_leaves           = TRUE,
  label_branch_points    = FALSE
)

# 按特定 meta 属性着色（例如 Genotype）
plot_cells(
  cds,
  color_cells_by        = "Genotype",
  show_trajectory_graph = TRUE
)

# --------------------------
# 9️⃣ 沿伪时间的基因表达变化
# --------------------------
AFD_genes <- rev(c("Tmem119", "Trem2", "Apoe", "Birc5", "Spp1", "Plac8"))

AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes, ]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)

plot_genes_in_pseudotime(
  AFD_lineage_cds,
  color_cells_by = "Stage",
  min_expr       = 0.5
)

# 也可以在 UMAP 上展示特定标记基因
plot_cells(
  cds,
  genes                  = rev(c("Tmem119", "Trem2", "Apoe", "Spp1", "Birc5", "Plac8")),
  label_cell_groups      = FALSE,
  show_trajectory_graph  = TRUE
)
```

---

If needed, you can split this file into `scRNAseq_QC_Integration.md` and `scRNAseq_Pseudotime.md`, but the current version keeps the full **single-cell pipeline** in one place for easier sharing and upload.
