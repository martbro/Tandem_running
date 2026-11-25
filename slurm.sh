#!/bin/bash
#SBATCH --job-name=acmg-snv
#SBATCH --cpus-per-task=1
#SBATCH --output=%x-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
exec "${ACMG_PYTHON:-python3}" -m acmg_snv vcf "$@"
