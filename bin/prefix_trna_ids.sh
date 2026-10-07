#!/usr/bin/env bash
# Prefix every ID= and Parent= of a tRNAscan-SE GFF3 with "<prefix>-", so that
# the tRNA IDs are unique next to the other features of the annotation
# (BRAKER4 rules/ncrna/run_trnascan.smk, commit 3535ed3).
# Copyright (c) 2025 Katharina Hoff. MIT License, see LICENSE-BRAKER4.
#
# Usage: prefix_trna_ids.sh <prefix> <tRNAscan-SE.gff3>  (writes to stdout)
set -euo pipefail
prefix=${1:?usage: prefix_trna_ids.sh <prefix> <file.gff3>}
gff=${2:?usage: prefix_trna_ids.sh <prefix> <file.gff3>}
awk -F'\t' -v OFS='\t' -v p="$prefix" '
    /^#/ { print; next }
    {
        $9 = ";" $9
        gsub(/;[[:space:]]*ID=/, ";ID=" p "-", $9)
        gsub(/;[[:space:]]*Parent=/, ";Parent=" p "-", $9)
        $9 = substr($9, 2)
        print
    }' "$gff"
