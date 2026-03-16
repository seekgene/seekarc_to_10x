#!/bin/bash

# ==============================================================================
# SeekARC_to_10X Pipeline Execution Script
# ==============================================================================
# This script demonstrates how to run the SeekARC_to_10X conversion pipeline
# for a Multiome (RNA+ATAC) sample. It sets up input file paths, output
# directories, and invokes the main python pipeline script.
#
# Usage:
#   bash run.sh
#
# Prerequisites:
#   1. Ensure the Conda environment is activated: conda activate seekarc_env
#   2. Verify input file paths exist.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Define Input Data Paths
# ------------------------------------------------------------------------------
# Replace these paths with your actual input files.
# RNA-seq Read 1 and Read 2
RNA_R1="./data/test_sample_GE_S180_L007_R1_001.fastq.gz"
RNA_R2="./data/test_sample_GE_S180_L007_R2_001.fastq.gz"

# ATAC-seq Read 1 and Read 2
ATAC_R1="./data/test_sample_arc_S33_L007_R1_001.fastq.gz"
ATAC_R2="./data/test_sample_arc_S33_L007_R2_001.fastq.gz"

# ------------------------------------------------------------------------------
# 2. Define Output Directory
# ------------------------------------------------------------------------------
# The directory where all results (QC, intermediate data, converted files) will be saved.
OUTDIR="./test_output"

# ------------------------------------------------------------------------------
# 3. Print Configuration (Optional)
# ------------------------------------------------------------------------------
echo "============================================================"
echo "Starting SeekARC_to_10X Pipeline Test"
echo "============================================================"
echo "Configuration:"
echo "  - RNA R1:    $RNA_R1"
echo "  - RNA R2:    $RNA_R2"
echo "  - ATAC R1:   $ATAC_R1"
echo "  - ATAC R2:   $ATAC_R2"
echo "  - Output Dir: $OUTDIR"
echo "============================================================"

# ------------------------------------------------------------------------------
# 4. Execute Pipeline
# ------------------------------------------------------------------------------
# Note: Ensure 'pipeline.py' is in the current directory or provide full path.
python pipeline.py \
    --mode "multiome" \
    --sample "test_sample" \
    --outdir "$OUTDIR" \
    --rna_r1 "$RNA_R1" \
    --rna_r2 "$RNA_R2" \
    --atac_r1 "$ATAC_R1" \
    --atac_r2 "$ATAC_R2" \
    --core 12

# ------------------------------------------------------------------------------
# 5. Completion
# ------------------------------------------------------------------------------
echo "============================================================"
echo "Pipeline execution finished."
echo "Please check results in: $OUTDIR"
echo "============================================================"