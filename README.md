# Han-Chen-Yang_Trem2_development

**A comprehensive bulk and single-cell RNA sequencing analysis pipeline for studying Trem2-dependent microglial development**  
**Tang Yi Laboratory, Xuanwu Hospital / Capital Medical University**

---

## 📌 Project Overview

This repository contains analysis documentation and pipelines associated with the project:

**“Han-Chen-Yang_Trem2_development”**  
developed in **Tang Yi Lab**, Department of Neurology, Xuanwu Hospital.

Microglia undergo dynamic developmental transitions across early postnatal stages (P7–P21).  
TREM2—a central microglial immune receptor—is essential not only in neurodegeneration, but also potentially in microglial maturation during development.

This repository summarizes the **bulk RNA-seq** and **single-cell RNA-seq** computational workflows used to analyze Trem2-dependent developmental trajectories.

---

## 📁 Repository Contents

```text
Han-Chen-Yang_Trem2_development/
│
├── bulkRNAseq_preprocess.md      # Raw FASTQ → gene counts workflow
├── bulkRNAseq1.md                # DEG, GO/KEGG, GSVA, ternary plot
├── bulkRNAseq2.md                # WGCNA, Mfuzz, UMAP-based signatures
│
└── scRNAseq_pipeline.md          # Single-cell QC, integration (scVI), pseudotime
```

> ⚠️ Note:  
> This repository only includes workflow documentation, not raw results or scripts.  
> Users can reference the Markdown files for step-by-step reproduction.

---

## 🧬 Bulk RNA-seq Workflow Summary

The bulk RNA-seq pipeline includes:

### Preprocessing

- MD5 checksum validation  
- FASTQ renaming & organization  
- Trim Galore QC  
- STAR genome alignment  
- RSEM quantification  
- Gene annotation & matrix merging  

### Downstream analyses

- Differential expression (DESeq2)  
- Developmental stage marker identification  
- Ternary composition analysis (P7–P14–P21)  
- GO / KEGG enrichment  
- GSVA pathway scoring  
- WGCNA co-expression networks  
- Mfuzz time-series clustering  

Detailed documentation:

- `bulkRNAseq_preprocess.md`  
- `bulkRNAseq1.md`  
- `bulkRNAseq2.md`  

---

## 🧫 Single-cell RNA-seq Workflow Summary

The single-cell pipeline covers:

### QC & Processing

- Per-sample Seurat object creation  
- DropletQC-derived metrics  
- Filtering based on expression, mt%, hb%, cell_probability, etc.  

### Integration

- Doublet removal  
- Batch-aware normalization  
- scVI-based latent embedding  
- UMAP visualization  
- Clustering and annotation  

### Trajectory Analysis

- Conversion to Monocle3  
- Graph learning  
- Pseudotime ordering  
- Lineage-specific gene dynamics  

Detailed documentation:

- `scRNAseq_pipeline.md`  

---

## 🧪 Biological Goal

This project aims to resolve:

- How Trem2 deficiency alters microglial transcriptional maturation  
- Whether developmental trajectories are delayed or altered  
- Which pathways and regulatory modules depend on TREM2  
- Stage-specific programs (P7 → P14 → P21) disrupted by Trem2-KO  
- Microglial lineage progression visualized at single-cell resolution  

---

## 🏛️ Laboratory

**Tang Yi Laboratory**  
Department of Neurology, Xuanwu Hospital  
Capital Medical University, Beijing, China  

Focus areas include:

- Microglial biology  
- Neuroimmune signaling  
- Alzheimer's disease & TREM2  
- Translational neuroimmunology  

---

## 👤 Maintainer

**Han-Chen Yang, M.D.**  
Department of Neurology  
Xuanwu Hospital, Capital Medical University  

Yusen Deng completed the bioinformatics analysis for this study.
