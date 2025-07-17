#!/usr/bin/env python3
import subprocess
import sys

jobid = sys.argv[1]
try:
    result = subprocess.check_output(
        ["sacct", "-j", jobid, "--format=State", "--noheader"],
        universal_newlines=True
    ).strip().split()[0]
except subprocess.CalledProcessError:
    result = "UNKNOWN"

if result.startswith("COMPLETED"):
    sys.exit(0)
elif result.startswith("FAILED") or result.startswith("CANCELLED") or result.startswith("TIMEOUT"):
    sys.exit(1)
else:
    sys.exit(2)

