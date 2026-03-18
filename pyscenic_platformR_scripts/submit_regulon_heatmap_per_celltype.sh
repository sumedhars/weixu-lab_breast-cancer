#!/bin/bash
#SBATCH --job-name=regulon_per_ct
#SBATCH --output=regulon_heatmap_logs/per_ct_%j.out
#SBATCH --error=regulon_heatmap_logs/per_ct_%j.err
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

# Column in adata.obs that has your cell type labels
CELL_TYPE_COL="cell_type"

# Column in adata.obs that distinguishes WT vs KO
CONDITION_COL="condition"

# Number of top regulons to select per cell type
N_TOP=15

# Output directory
OUTPUT_DIR="regulon_heatmap_per_celltype"

# =============================================================================
# RUN
# =============================================================================

echo "============================================"
echo "Per-Cell-Type Regulon Heatmap: WT vs KO"
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
echo "  Per cell type outputs (in <OUTPUT_DIR>/<celltype>/):"
echo "    - condition_tvalues.csv              (all regulon t-values)"
echo "    - regulon_tvalue_heatmap.png         (top N condition t-values)"
echo "    - regulon_wt_vs_ko_heatmap.png       (mean AUCell: WT vs KO rows)"
echo "    - regulon_violin_aucell.png          (AUCell distributions split by WT/KO)"
echo "    - regulon_dotplot_tf_expression.png  (TF gene expression dotplot by WT/KO)"
echo ""
echo "  Summary outputs (in <OUTPUT_DIR>/):"
echo "    - condition_tvalues_all_celltypes.csv (combined t-value matrix)"
echo "    - summary_condition_heatmap.png       (all cell types, condition effect)"
echo ""

${PYTHON_BIN} plot_regulon_heatmap_per_celltype.py \
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
