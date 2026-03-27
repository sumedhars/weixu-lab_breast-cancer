#!/bin/bash
#SBATCH --job-name=scenic_diag
#SBATCH --output=pyscenic_logs/diagnose_%j.out
#SBATCH --error=pyscenic_logs/diagnose_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p pyscenic_logs

# =============================================================================
# ENVIRONMENT
# =============================================================================

eval "$(conda shell.bash hook)"
conda activate pyscenic
export PYTHONNOUSERSITE=1

PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/pyscenic/bin/python"

# =============================================================================
# PATHS — UPDATE THESE TO MATCH YOUR OUTPUT DIRECTORY
# =============================================================================

# Path to modules.pkl from your completed Phase I run
MODULES_PKL="regulon_heatmap_output_all/modules.pkl"

# Path to adjacencies.tsv (optional but recommended)
ADJACENCIES_TSV="regulon_heatmap_output_all/adjacencies.tsv"

# Same resource paths as your main pipeline
DB_DIR="pyscenic_resources/databases"
MOTIF_ANNOTATIONS="pyscenic_resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# =============================================================================
# RUN DIAGNOSTIC
# =============================================================================

echo "============================================"
echo "SCENIC TF Regulon Diagnostic"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start:  $(date)"
echo "============================================"

${PYTHON_BIN} diagnose_tf_modules.py \
    --modules_pkl "${MODULES_PKL}" \
    --db_dir "${DB_DIR}" \
    --motif_annotations "${MOTIF_ANNOTATIONS}" \
    --tfs Tfap2b Foxa1 Esr1 \
    --adjacencies_tsv "${ADJACENCIES_TSV}"

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
