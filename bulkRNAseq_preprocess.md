# Tangyilab Bulk RNA-seq Preprocessing Pipeline

This document summarizes the **raw data → gene-level count matrix** preprocessing steps for bulk RNA-seq, including:

1. Master pipeline launcher (`autoRNA.sh`)  
2. Step 1.1–1.2: MD5 check and FASTQ collection / renaming (`CheckandCollect.sh`)  
3. Step 2: QC and adapter trimming with Trim Galore (`QC.sh`)  
4. Step 3: Genome alignment with STAR (`STAR_Mapping.sh`)  
5. Step 4: Gene-level quantification with RSEM (`RSEM_Count.sh`)  
6. Step 5: Merging RSEM `expected_count` matrices and gene annotation (`merge_rsem_expected_count_with_annot.py`)

---

## 1. Master Bulk RNA-seq Pipeline Launcher (`autoRNA.sh`)

```bash
#!/bin/bash

echo "#############################################"
echo "###                                       ###"
echo "###   #Tangyilab_auto_bulk_pipline V1.0   ###"
echo "###                                       ###"
echo "#############################################"


echo "#############################################"
echo "###                                       ###"
echo "###  created by Deng YuSen and ZhouJing   ###"
echo "###                                       ###"
echo "#############################################"

CALL_DIR=$(pwd)


bash /data_result/dengys/bulk/autobulk/CheckandCollect.sh "$CALL_DIR" 

echo "#############################################"
echo "###                                       ###"
echo "###          step 2 : QC                  ###"
echo "###                                       ###"
echo "#############################################"
bash /data_result/dengys/bulk/autobulk/QC.sh "$CALL_DIR" > "$CALL_DIR/QC.log"

echo "#############################################"
echo "###                                       ###"
echo "###          QC finished                  ###"
echo "###                                       ###"
echo "#############################################"

echo "#############################################"
echo "###                                       ###"
echo "###          step 3 : MApping             ###"
echo "###                                       ###"
echo "#############################################"

bash /data_result/dengys/bulk/autobulk/STAR_Mapping.sh "$CALL_DIR" > "$CALL_DIR/Mapping.log"

echo "#############################################"
echo "###                                       ###"
echo "###          MApping finished             ###"
echo "###                                       ###"
echo "#############################################"


echo "#############################################"
echo "###                                       ###"
echo "###          step 4 : Counting            ###"
echo "###                                       ###"
echo "#############################################"

bash /data_result/dengys/bulk/autobulk/RSEM_Count.sh "$CALL_DIR" > "$CALL_DIR/Counting.log"

echo "#############################################"
echo "###                                       ###"
echo "###          Counting finished            ###"
echo "###                                       ###"
echo "#############################################"


echo "#############################################"
echo "###                                       ###"
echo "###          step 5 : Merging             ###"
echo "###                                       ###"
echo "#############################################"


/home/dengys/anaconda3/bin/python /data_result/dengys/bulk/autobulk/merge_rsem_expected_count_with_annot.py \
 --pattern "$CALL_DIR/RSEM_quant/*.genes.results" \
 --annot "/data_result/dengys/bulk/autobulk/ref/extracted_gene_info_mmu.csv" \
 --out-prefix "expected_count" \
 --out-dir "$CALL_DIR"

echo "#############################################"
echo "###                                       ###"
echo "###          Merging finished             ###"
echo "###                                       ###"
echo "#############################################"
```

---

## 2. Step 1.1–1.2: MD5 Check & Collect FASTQ Files (`CheckandCollect.sh`)

```bash
#!/bin/bash

# Stylized header
echo "#############################################"
echo "###                                       ###"
echo "###             Step:1.1                  ###"
echo "###   MD5 Check and Archive Script        ###"
echo "###                                       ###"
echo "#############################################"

LOG_FILE="$1/md5.log"

# 如果 log 存在
if [ -f "$LOG_FILE" ]; then
    last_line=$(tail -n 1 "$LOG_FILE")
    if [[ "$last_line" == "all file is right" ]]; then
        echo "MD5 log already verified. Skipping MD5 check."
        VERIFY_OK=true
    else
        echo "File break check and redownload"
        exit 1
    fi
else
    # log 不存在，执行 MD5 校验
    md5_files=( "$1"/*.md5 )
    if [ -e "${md5_files[0]}" ]; then
        echo "Found .md5 files in $1. Starting verification..."
        > "$LOG_FILE"  # 清空/新建日志
        for md5_file in "${md5_files[@]}"; do
            echo "Verifying using $md5_file"
            md5sum -c "$md5_file" >> "$LOG_FILE" 2>&1
        done

        # 检查是否所有文件都通过
        if grep -q "FAILED" "$LOG_FILE"; then
            echo "File break check and redownload"
            exit 1
        else
            echo "all file is right" >> "$LOG_FILE"
            echo "MD5 verification completed successfully."
            VERIFY_OK=true
        fi
    else
        echo "No .md5 files found in $1."
        exit 1
    fi
fi

echo "Step 1 completed."

#################################################
# Step 2: Collect all files into a tar.gz archive
#################################################

#!/usr/bin/env bash
set -euo pipefail

if [ "$VERIFY_OK" = true ]; then
    echo "Start file organizing..."
else
    exit 1
fi

echo "#############################################"
echo "###                                       ###"
echo "###             Step:1.2                  ###"
echo "###         Collect all files             ###"
echo "###                                       ###"
echo "#############################################"

LOG_FILE_COL="$1/collect.log"
RAW_DIR="$1/raw"
mkdir -p "$RAW_DIR"

# 如果 log 文件存在并且最后一行是 ALL_DONE，说明已经跑过
if [ -f "$LOG_FILE_COL" ]; then
    last_line=$(tail -n 1 "$LOG_FILE_COL")
    if [[ "$last_line" == "ALL_DONE" ]]; then
        echo "Fastq file is already collected. Skipping."
        exit 0
    fi
fi

# 统一命名规则
normalize_name() {
    local newname="$1"
    newname=$(echo "$newname" | sed -E 's/\.R1\.fq\.gz$/_R1.fq.gz/')
    newname=$(echo "$newname" | sed -E 's/\.R2\.fq\.gz$/_R2.fq.gz/')
    newname=$(echo "$newname" | sed -E 's/R_1\.fq\.gz$/_R1.fq.gz/')
    newname=$(echo "$newname" | sed -E 's/R_2\.fq\.gz$/_R2.fq.gz/')
    newname=$(echo "$newname" | sed -E 's/_1\.fq\.gz$/_R1.fq.gz/')
    newname=$(echo "$newname" | sed -E 's/_2\.fq\.gz$/_R2.fq.gz/')
    echo "$newname"
}

# 逐个文件 copy+校验+删除+写 log
process_file() {
    local f="$1"
    local base=$(basename "$f")
    local newname
    newname=$(normalize_name "$base")
    local dst="$RAW_DIR/$newname"
    local tmp="$dst.tmp"

    # 如果 log 已记录该文件 OK，跳过
    if grep -Fq "OK $base -> $newname" "$LOG_FILE_COL" 2>/dev/null; then
        echo "Skip (already logged): $f"
        return
    fi

    echo "Copying $f -> $dst"
    cp "$f" "$tmp"
    if cmp -s "$f" "$tmp"; then
        mv "$tmp" "$dst"
        rm -f "$f"
        echo "OK $base -> $newname $(stat -c %s "$dst") $(date +%s)" >> "$LOG_FILE_COL"
    else
        echo "ERROR: Copy verify failed for $f"
        rm -f "$tmp"
        exit 1
    fi
}

echo "Start file collecting..."
shopt -s nullglob
for f in "$1"/*.fq.gz; do
    [ -e "$f" ] || continue
    process_file "$f"
done

echo "ALL_DONE" >> "$LOG_FILE_COL"
echo "File organization completed. All FASTQ files copied to $RAW_DIR and verified."
```

---

## 3. Step 2: QC and Adapter Trimming with Trim Galore (`QC.sh`)

```bash
#!/bin/bash
# Stylized header
echo "#############################################"
echo "###                                       ###"
echo "###              Step: 2                  ###"
echo "###         QC with Trim Galore           ###"
echo "###                                       ###"
echo "#############################################"

# Set directories
RAW_DIR="$1/raw"
CLEAN_DIR="$1/clean"
mkdir -p "$CLEAN_DIR"

# Loop through R1 files
for r1 in "$RAW_DIR"/*R1.fq.gz; do
    r2="${r1/R1.fq.gz/R2.fq.gz}"
    prefix=$(basename "$r1" R1.fq.gz)

    # Check if output files already exist
    if [ -f "$CLEAN_DIR/${prefix}_R1_val_1.fq.gz" ] && [ -f "$CLEAN_DIR/${prefix}_R2_val_2.fq.gz" ]; then
        echo "⏭️ Skipping $prefix (already trimmed)"
        continue
    fi

    echo "▶️ Trimming $prefix ..."
    /data_result/dengys/bulk/autobulk/software/TrimGalore-0.6.10/trim_galore --paired \
        --gzip \
        --path_to_cutadapt /home/dengys/anaconda3/envs/bulk/bin/cutadapt \
        --cores 40 \
        --fastqc \
        --output_dir "$CLEAN_DIR" \
        "$r1" "$r2"
    echo "✅ Done trimming $prefix"
done
```

---

## 4. Step 3: Genome Alignment with STAR (`STAR_Mapping.sh`)

```bash
#!/bin/bash

echo "#############################################"
echo "###                                       ###"
echo "###              Step: 3                  ###"
echo "###            STAR mapping               ###"
echo "###                                       ###"
echo "#############################################"
ulimit -n 65535
# 参数说明
CALL_DIR=$1  # 输入目录，包含FASTQ文件
# 自动检测核心数并设置线程数
TOTAL_CORES=$(nproc)
THREADS=$(( TOTAL_CORES / 2 ))
if [ "$THREADS" -gt 40 ]; then
    THREADS=40
fi
echo "🧠 Detected $TOTAL_CORES cores, using $THREADS threads."
GENOME_DIR="/data_result/dengys/bulk/autobulk/ref/STARindex"  # STAR索引路径
OUT_DIR="${CALL_DIR}/STAR_output"  # 输出目录
# 创建输出目录
mkdir -p "$OUT_DIR"

# 遍历FASTQ文件进行比对
for fq1 in ${CALL_DIR}/clean/*_R1_val_1.fq.gz; do
    fq2=${fq1/_R1_val_1.fq.gz/_R2_val_2.fq.gz}
    sample=$(basename "$fq1" _R1_val_1.fq.gz)

    echo "▶️ Mapping sample: $sample"

    /home/dengys/anaconda3/envs/bulk/bin/STAR \
        --runThreadN $THREADS \
        --genomeDir $GENOME_DIR \
        --readFilesIn $fq1 $fq2 \
        --readFilesCommand zcat \
        --outFileNamePrefix ${OUT_DIR}/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts
    echo "✅ Finished mapping: $sample"
done

echo "🎉 All samples mapped successfully!"
```

---

## 5. Step 4: Gene-level Quantification with RSEM (`RSEM_Count.sh`)

```bash
#!/bin/bash

echo "#############################################"
echo "###                                       ###"
echo "###              Step: 4                  ###"
echo "###            RSEM mapping               ###"
echo "###                                       ###"
echo "#############################################"
ulimit -n 65535
# 参数说明
CALL_DIR=$1  # 输入目录，包含FASTQ文件
# 自动检测核心数并设置线程数
TOTAL_CORES=$(nproc)
THREADS=$(( TOTAL_CORES / 2 ))
if [ "$THREADS" -gt 40 ]; then
    THREADS=40
fi
echo "🧠 Detected $TOTAL_CORES cores, using $THREADS threads."


# 路径配置
RSEM_INDEX="/data_result/dengys/bulk/autobulk/ref/RESMindex/resmindex"
INPUT_DIR="${CALL_DIR}/STAR_output" 
OUTPUT_DIR="${CALL_DIR}/RSEM_quant"
mkdir -p "$OUTPUT_DIR"



# 遍历所有 _Aligned.toTranscriptome.out.bam 文件
find "$INPUT_DIR" -type f -name "*_Aligned.toTranscriptome.out.bam" | while read -r bam; do
    filename=$(basename "$bam")
    prefix="${filename%%_Aligned.toTranscriptome.out.bam}"

    echo "🧮 Running RSEM quantification for $prefix..."

    /home/dengys/anaconda3/envs/bulk/bin/rsem-calculate-expression \
        --alignments \
        --bam \
        --paired-end \
        --no-bam-output \
        --estimate-rspd \
        --strandedness none \
        --num-threads "$THREADS" \
        "$bam" \
        "$RSEM_INDEX" \
        "$OUTPUT_DIR/$prefix"

    echo "✅ Done $prefix"
done

echo "🎉 All samples counted successfully!"
```

---

## 6. Step 5: Merge RSEM `expected_count` Matrices and Annotate Genes (`merge_rsem_expected_count_with_annot.py`)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge RSEM *.genes.results expected_count, annotate with gene info, filter to protein_coding,
and collapse multiple ENS IDs mapping to the same Gene_Name.

新增:
- --out-dir  指定输出目录；默认=检测到的 RSEM_quant 上一级目录。

Usage 示例：
  python merge_rsem_expected_count_with_annot.py \
    --pattern "*.genes.results" \
    --annot "/data_result/dengys/bulk/autobulk/ref/extracted_gene_info_mmu.csv" \
    --out-prefix "expected_count" \
    --out-dir "/data_result/dengys/bulk/2509_YHC_microgliadevelopment/YHC_Bulk_Microglia_develop"

Author: dengys
Date: 2025
Description: Merges RSEM gene results across samples and performs gene annotation filtering

Input: 
  - *.genes.results files from RSEM quantification
  - Gene annotation CSV file

Output:
  - Raw merged matrix
  - Filtered and annotated matrix
  - Gene name collapsed matrix

Dependencies: pandas
"""

import argparse
import glob
import os
import sys
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def infer_sample_name(path: str) -> str:
    """
    Infer sample name from file path.

    Parameters
    ----------
    path : str
        Path to genes.results file

    Returns
    -------
    str
        Sample name derived from filename
    """
    base = os.path.basename(path)
    if base.endswith(".genes.results"):
        return base[:-len(".genes.results")]
    else:
        return os.path.splitext(base)[0]


def read_one_genes_results(file_path: str, value_col: str = "expected_count") -> pd.DataFrame:
    """
    Read single RSEM genes.results file.

    Parameters
    ----------
    file_path : str
        Path to genes.results file
    value_col : str, default "expected_count"
        Column to extract values from

    Returns
    -------
    pd.DataFrame
        DataFrame with gene_id and value columns

    Raises
    ------
    ValueError
        If required columns are missing from file
    """
    try:
        df = pd.read_csv(file_path, sep="\t", engine="python", on_bad_lines="skip")
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        raise

    # Validate required columns
    required = ["gene_id", "transcript_id(s)", value_col]
    missing_cols = [c for c in required if c not in df.columns]
    if missing_cols:
        raise ValueError(f"{file_path} 缺少列: {missing_cols}")

    # Extract and clean data
    out = df[["gene_id", value_col]].copy()
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")

    # Check for conversion issues
    if out[value_col].isna().sum() > 0:
        logger.warning(f"Found {out[value_col].isna().sum()} non-numeric values in {file_path}")

    return out


def merge_all_samples(pattern: str, value_col: str = "expected_count") -> tuple[pd.DataFrame, list]:
    """
    Merge all RSEM results files matching the pattern.

    Parameters
    ----------
    pattern : str
        Glob pattern to match files
    value_col : str, default "expected_count"
        Column to extract from each file

    Returns
    -------
    tuple[pd.DataFrame, list]
        Merged dataframe and list of processed files

    Raises
    ------
    SystemExit
        If no files match the pattern
    """
    files = sorted(glob.glob(pattern))
    if not files:
        logger.error(f"未找到匹配文件：{pattern}")
        sys.exit(1)

    logger.info(f"Found {len(files)} files matching pattern: {pattern}")

    merged = None
    for fp in files:
        sample_name = infer_sample_name(fp)
        logger.info(f"Processing {fp} -> {sample_name}")

        df = read_one_genes_results(fp, value_col)
        df = df.rename(columns={value_col: sample_name})

        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on="gene_id", how="outer")

    # Fill missing values and set index
    merged = merged.fillna(0.0)
    merged = merged.set_index("gene_id").sort_index()

    logger.info(f"Merged matrix shape: {merged.shape}")
    return merged, files


def detect_default_output_dir(files: list) -> str:
    """
    Detect default output directory from file paths.

    If all files have common path ending with 'RSEM_quant', 
    return the parent directory, otherwise return common path.

    Parameters
    ----------
    files : list
        List of file paths

    Returns
    -------
    str
        Default output directory path
    """
    common_base = os.path.commonpath([os.path.abspath(f) for f in files])

    # Handle case where common_base is a file
    if os.path.isfile(common_base):
        common_base = os.path.dirname(common_base)

    # Check if last directory component is 'RSEM_quant'
    tail = os.path.basename(common_base.rstrip(os.sep))
    if tail == "RSEM_quant":
        return os.path.dirname(common_base)
    else:
        return common_base


def load_gene_annotation(annot_csv: str) -> pd.DataFrame:
    """
    Load and validate gene annotation file.

    Parameters
    ----------
    annot_csv : str
        Path to annotation CSV file

    Returns
    -------
    pd.DataFrame
        Gene annotation dataframe

    Raises
    ------
    ValueError
        If required columns are missing
    FileNotFoundError
        If annotation file doesn't exist
    """
    if not os.path.exists(annot_csv):
        raise FileNotFoundError(f"Annotation file not found: {annot_csv}")

    try:
        annot = pd.read_csv(annot_csv)
    except Exception as e:
        logger.error(f"Error reading annotation file {annot_csv}: {e}")
        raise

    # Validate required columns
    required_cols = ["Gene_ID", "Gene_Type", "Gene_Name"]
    missing_cols = [c for c in required_cols if c not in annot.columns]
    if missing_cols:
        raise ValueError(f"注释文件缺少列: {missing_cols}")

    # Remove duplicates based on Gene_ID
    initial_rows = len(annot)
    annot = annot.drop_duplicates(subset=["Gene_ID"])
    final_rows = len(annot)

    if initial_rows != final_rows:
        logger.warning(f"Removed {initial_rows - final_rows} duplicate Gene_IDs from annotation")

    logger.info(f"Loaded annotation for {len(annot)} genes")
    return annot


def create_output_path(out_dir: str, filename: str) -> str:
    """Create full output file path."""
    return os.path.join(out_dir, filename)


def main():
    """Main analysis workflow."""
    # Parse arguments
    parser = argparse.ArgumentParser(description="Merge RSEM results with gene annotation")
    parser.add_argument("--pattern", default="*.genes.results", 
                       help="RSEM genes.results 匹配模式")
    parser.add_argument("--value-col", default="expected_count", 
                       choices=["expected_count", "TPM", "FPKM"],
                       help="汇总列（默认 expected_count）")
    parser.add_argument("--annot", required=True, 
                       help="基因注释 CSV（需含 Gene_ID,Gene_Type,Gene_Name）")
    parser.add_argument("--out-prefix", default="expected_count", 
                       help="输出文件名前缀")
    parser.add_argument("--out-dir", default=None, 
                       help="输出目录（默认=RSEM_quant 上一级目录）")
    parser.add_argument("--keep-non-protein", action="store_true", 
                       help="不过滤非蛋白编码（默认会过滤）")
    parser.add_argument("--no-collapse", action="store_true", 
                       help="不按 Gene_Name 聚合（默认会聚合）")
    parser.add_argument("--agg", default="sum", 
                       choices=["sum", "mean", "median", "max"],
                       help="按 Gene_Name 聚合方式（默认 sum）")

    args = parser.parse_args()

    logger.info("Starting RSEM results merging and annotation")

    try:
        # Step 1: Merge all samples
        logger.info("Step 1: Merging all sample files")
        raw_matrix, files = merge_all_samples(args.pattern, args.value_col)

        # Step 2: Determine output directory
        out_dir = args.out_dir or detect_default_output_dir(files)
        os.makedirs(out_dir, exist_ok=True)
        logger.info(f"Output directory: {out_dir}")

        # Step 3: Save raw matrix
        raw_output_path = create_output_path(out_dir, f"{args.out_prefix}.raw_matrix.csv")
        raw_matrix.to_csv(raw_output_path)
        logger.info(f"[输出] 原始合并矩阵：{raw_output_path}  形状={raw_matrix.shape}")

        # Step 4: Load annotation and merge
        logger.info("Step 4: Loading gene annotation and filtering")
        annot = load_gene_annotation(args.annot)

        # Prepare matrix for annotation
        mat = raw_matrix.reset_index().rename(columns={"index": "gene_id"})
        mat = mat.merge(annot, left_on="gene_id", right_on="Gene_ID", how="left")

        # Log annotation statistics
        annotated_genes = mat["Gene_Name"].notna().sum()
        total_genes = len(mat)
        logger.info(f"Successfully annotated {annotated_genes}/{total_genes} genes")

        # Step 5: Filter to protein coding if requested
        if not args.keep_non_protein:
            before_filter = mat.shape[0]
            mat = mat[mat["Gene_Type"] == "protein_coding"].copy()
            after_filter = mat.shape[0]
            logger.info(f"[过滤] 保留 protein_coding：{before_filter} -> {after_filter}")

        # Prepare gene ID level output
        sample_cols = [c for c in mat.columns 
                      if c not in {"gene_id", "Gene_ID", "Gene_Type", "Gene_Name"}]

        id_level = mat[["gene_id", "Gene_Name", "Gene_Type"] + sample_cols].copy()
        id_level = id_level.set_index("gene_id").sort_index()

        # Save gene ID level results
        id_filtered_path = create_output_path(out_dir, f"{args.out_prefix}.geneid_annot.filtered.csv")
        id_level.to_csv(id_filtered_path)
        logger.info(f"[输出] 基因ID层面（已注释/过滤）：{id_filtered_path}  形状={id_level.shape}")

        # Step 6: Gene name aggregation
        if not args.no_collapse:
            logger.info(f"Step 6: Collapsing by Gene_Name using {args.agg} aggregation")

            # Prepare data for aggregation
            collapse_data = mat[["Gene_Name"] + sample_cols].copy()
            collapse_data = collapse_data.dropna(subset=["Gene_Name"])

            if len(collapse_data) == 0:
                logger.warning("No genes with valid Gene_Name found for collapsing")
            else:
                # Perform aggregation
                collapsed = (collapse_data
                           .groupby("Gene_Name", as_index=True)
                           .agg(args.agg)
                           .sort_index())

                # Save collapsed results
                collapsed_path = create_output_path(out_dir, 
                                                  f"{args.out_prefix}.genename_collapsed.filtered.csv")
                collapsed.to_csv(collapsed_path)
                logger.info(f"[输出] 基因名聚合矩阵（{args.agg}）：{collapsed_path}  形状={collapsed.shape}")

        logger.info("[完成] Analysis completed successfully")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
```
