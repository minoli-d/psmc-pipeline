#!/bin/bash

ASSEMBLY_ROOT="/global/scratch/users/nicolas931010/sv_diversity/assemblies"

for accession_dir in "$ASSEMBLY_ROOT"/*; do
    accession=$(basename "$accession_dir")
    jsonl="$accession_dir/sequence_report.jsonl"

    if [[ -f "$jsonl" ]]; then
        # Collect all chrNames for chromosomes, remove those containing X, Y, Z, or W (case-insensitive)
        autosomes=$(jq -R 'fromjson? | select(. != null) |
            select(.assignedMoleculeLocationType == "Chromosome") | .chrName' "$jsonl" \
            | tr -d '"' \
            | grep -Ev '[XYZWxyzw]' \
            | sort -u \
            | paste -sd "," -)
        echo -e "${accession}\t${autosomes}"
    else
        echo -e "${accession}\tNA"
    fi
done > autosomes_by_accession.tsv

