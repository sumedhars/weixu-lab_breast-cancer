#!/bin/bash
#SBATCH --job-name=pyscenic_run
#SBATCH --output=pyscenic_logs/pyscenic_%A_%a.out
#SBATCH --error=pyscenic_logs/pyscenic_%A_%a.err
#SBATCH --time=128:00:00
#SBATCH --mem=950G
#SBATCH --cpus-per-task=512
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p pyscenic_logs

# =============================================================================
# ENVIRONMENT — Isolate conda env from ~/.local site-packages
# =============================================================================

eval "$(conda shell.bash hook)"
conda activate pyscenic

# CRITICAL: Prevent ~/.local/lib/python*/site-packages from leaking into
# the conda environment. Without this, stale or conflicting packages
# installed via `pip install --user` will shadow conda-managed packages.
export PYTHONNOUSERSITE=1

PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/pyscenic/bin/python"

# =============================================================================
# PATHS — UPDATE THESE TO MATCH YOUR FILE LOCATIONS
# =============================================================================

# Input h5ad file (combined WT + KO)
H5AD="/mnt/scratch/group/hqdinh2/weixu_lab_pyscenic/CTR9_snRNASeq/CTR9_snRNASeq_full.h5ad"

# Column in adata.obs that distinguishes conditions, and the labels
CONDITION_COL="sample"
WT_LABEL="WT_DM"
KO_LABEL="KO_DM"

TF_LIST="pyscenic_resources/allTFs_mm.txt"
DB_DIR="pyscenic_resources/databases"
MOTIF_ANNOTATIONS="pyscenic_resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# Output
OUTPUT_DIR="pyscenic_output"

# Number of workers for GRNBoost2 and AUCell (set to match --cpus-per-task or lower)
N_WORKERS=512

# =============================================================================
# RUN PIPELINE
# =============================================================================

echo "============================================"
echo "pySCENIC Pipeline — CTR9 KO vs WT"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "PYTHONNOUSERSITE: ${PYTHONNOUSERSITE}"
echo "Python: ${PYTHON_BIN}"
echo "Start:  $(date)"
echo "============================================"

${PYTHON_BIN} run_pyscenic.py \
    --h5ad "${H5AD}" \
    --condition_col "${CONDITION_COL}" \
    --wt_label "${WT_LABEL}" \
    --ko_label "${KO_LABEL}" \
    --tf_list "${TF_LIST}" \
    --db_dir "${DB_DIR}" \
    --motif_annotations "${MOTIF_ANNOTATIONS}" \
    --output_dir "${OUTPUT_DIR}" \
    --n_workers ${N_WORKERS}

# Uncomment below to resume from a checkpoint if a previous run partially completed:
# ${PYTHON_BIN} run_pyscenic.py \
#     --h5ad "${H5AD}" \
#     --condition_col "${CONDITION_COL}" \
#     --wt_label "${WT_LABEL}" \
#     --ko_label "${KO_LABEL}" \
#     --tf_list "${TF_LIST}" \
#     --db_dir "${DB_DIR}" \
#     --motif_annotations "${MOTIF_ANNOTATIONS}" \
#     --output_dir "${OUTPUT_DIR}" \
#     --n_workers ${N_WORKERS} \
#     --resume_from regulons

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
