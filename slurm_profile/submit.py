#!/usr/bin/env python

import os
import sys
import subprocess
from snakemake.utils import read_job_properties

# Load job properties
jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# Extract resource info
rule = job_properties["rule"]
threads = job_properties.get("threads", 1)
resources = job_properties.get("resources", {})
mem = int(resources.get("mem_mb", 4000))
time = resources.get("time", "01:00:00")

# Wildcard string for logging
wildcards = job_properties.get("wildcards", {})
wildcard_str = ".".join(f"{k}-{v}" for k, v in wildcards.items()) if wildcards else "no_wc"

# Make log directory
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)
log_out = os.path.join(log_dir, f"{rule}.{wildcard_str}.out")
log_err = os.path.join(log_dir, f"{rule}.{wildcard_str}.err")

# Build sbatch command
sbatch_cmd = [
    "sbatch",
    f"--job-name=smk-{rule}",
    f"--output={log_out}",
    f"--error={log_err}",
    f"--time={time}",
    f"--mem={mem}",
    f"--cpus-per-task={threads}",
    "--nodes=1",
    "--ntasks=1",
    "--partition=savio4_htc",
    "--qos=savio_lowprio",
    "--account=co_genomicdata",
    "--exclude=n0074.savio4",
    jobscript
]

# Submit the job
subprocess.run(sbatch_cmd)

