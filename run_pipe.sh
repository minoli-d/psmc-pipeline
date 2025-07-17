#!/bin/bash

snakemake \
  --profile slurm_profile \
  --executor slurm \
  --jobs 100
