#!/bin/bash

set -euo pipefail

echo "Starting job $SLURM_JOB_ID on host $(hostname)"
echo "Working directory: $(pwd)"

source activate psmc_env

{exec_job}

