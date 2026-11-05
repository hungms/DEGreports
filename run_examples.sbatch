#!/bin/bash
# Run the two example pipelines (one array task each) with the project renv.
#   sbatch run_examples.sbatch          # both
#   sbatch --array=1 run_examples.sbatch  # just one (0-1, see SCRIPTS below)
# Outputs: example_run_<dataset>_<engine>_<mode>/ ; console log in its logs/ folder.
#SBATCH --job-name=degreports-examples
#SBATCH --partition=ncpu
#SBATCH --ntasks=30
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-1
#SBATCH --chdir=/nemo/lab/caladod/working/Matthew/git-repos/r-dev/DEGreports
#SBATCH --output=example_run_logs/%x_%A_%a.out

set -euo pipefail

SCRIPTS=(run_gse189410_deseq2_wald.R run_gse189410_deseq2_lrt.R)
OUTDIRS=(example_run_GSE189410_deseq2_wald example_run_GSE189410_deseq2_lrt)
script=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
outdir=${OUTDIRS[$SLURM_ARRAY_TASK_ID]}

# pandoc for the HTML report (not in the seurat5 env).
module load Pandoc/2.19

# R 4.3.2 from the seurat5 conda env; packages come from the project renv
# (activated by .Rprofile -> renv/activate.R; .Renviron sets the renv cache).
RSCRIPT=/nemo/lab/caladod/working/Matthew/.conda/envs/seurat5/bin/Rscript

# Keep BLAS/OpenMP single-threaded: parallelism comes from render workers
# (= SLURM_NTASKS), and nested threads per forked worker oversubscribe the node.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p "$outdir/logs"
echo "[$(date)] $script on $(hostname) with $SLURM_NTASKS tasks -> $outdir"
"$RSCRIPT" "$script" > "$outdir/logs/console.log" 2>&1
echo "[$(date)] done: $outdir/index.html"
