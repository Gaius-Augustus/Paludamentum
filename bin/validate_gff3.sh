#!/usr/bin/env bash
# Validate a published GFF3 (contract item 6, docs/postprocessing.md):
# GenomeTools gt gff3validator must accept it, and gffread must print no
# warning about missing parents, overlapping exons or invalid CDS.
# Exit 1 with the validator output on failure.
#
# Usage: validate_gff3.sh <file.gff3>
set -uo pipefail

gff=${1:?usage: validate_gff3.sh <file.gff3>}
log=$(mktemp)
trap 'rm -f "$log"' EXIT

if ! gt gff3validator "$gff" > "$log" 2>&1; then
    echo "validate_gff3.sh: gt gff3validator rejects $gff:" >&2
    head -n 50 "$log" >&2
    exit 1
fi

# -E: expose (warn about) duplicate IDs, missing parents, invalid CDS
gffread -E "$gff" -o /dev/null > "$log" 2>&1
if grep -i -E 'warning|error|discarded|invalid|overlap|no parent|not found' "$log" > /dev/null; then
    echo "validate_gff3.sh: gffread -E warns about $gff:" >&2
    head -n 50 "$log" >&2
    exit 1
fi
echo "validate_gff3.sh: $gff is valid"
