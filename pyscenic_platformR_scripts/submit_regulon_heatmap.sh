#!/bin/bash
#SBATCH --job-name=regulon_heatmap
#SBATCH --output=regulon_heatmap_logs/heatmap_%j.out
#SBATCH --error=regulon_heatmap_logs/heatmap_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p regulon_heatmap_logs

# =============================================================================
# ENVIRONMENT
# =============================================================================

eval "$(conda shell.bash hook)"
conda activate pyscenic

export PYTHONNOUSERSITE=1

PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/pyscenic/bin/python"

# =============================================================================
# PATHS — UPDATE THESE TO MATCH YOUR FILE LOCATIONS
# =============================================================================

# The h5ad output from your pySCENIC run (contains adata.obs with cell types)
H5AD="pyscenic_output/adata_with_aucell.h5ad"

# The AUCell matrix CSV from pySCENIC (cells x regulons)
AUC_CSV="pyscenic_output/auc_matrix.csv"

# Column in adata.obs that has your 8 cell type labels
CELL_TYPE_COL="cell_type"

# Column in adata.obs that distinguishes WT vs KO (covariate in the linear model)
CONDITION_COL="condition"

# Number of top regulons to select per cell type
N_TOP=5

# Output directory
OUTPUT_DIR="regulon_heatmap_output"

# =============================================================================
# RUN
# =============================================================================

echo "============================================"
echo "SCENIC Regulon Heatmap — Figure 2D Style"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "Start:  $(date)"
echo "============================================"
echo ""
echo "  H5AD:           ${H5AD}"
echo "  AUC_CSV:        ${AUC_CSV}"
echo "  CELL_TYPE_COL:  ${CELL_TYPE_COL}"
echo "  CONDITION_COL:  ${CONDITION_COL}"
echo "  N_TOP:          ${N_TOP}"
echo "  OUTPUT_DIR:     ${OUTPUT_DIR}"
echo ""
echo "  Outputs:"
echo "    - regulon_heatmap_clustered.png        (original clustered heatmap)"
echo "    - regulon_heatmap_sorted.png           (sorted by top cell type)"
echo "    - regulon_heatmap_condition_sorted.png  (celltype_condition rows, sorted)"
echo "    - regulon_heatmap_condition_clustered.png (celltype_condition rows, clustered)"
echo ""

${PYTHON_BIN} plot_regulon_heatmap.py \
    --h5ad "${H5AD}" \
    --auc_csv "${AUC_CSV}" \
    --cell_type_col "${CELL_TYPE_COL}" \
    --condition_col "${CONDITION_COL}" \
    --n_top ${N_TOP} \
    --output_dir "${OUTPUT_DIR}"

echo ""
echo "============================================"
echo "Finished: $(date)"
echo "============================================"
