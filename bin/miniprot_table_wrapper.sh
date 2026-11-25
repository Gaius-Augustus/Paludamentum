#!/usr/bin/env bash
# compleasm and BUSCO in genome mode find their genes with miniprot and take
# the protein sequences from its --trans output, but have no option for the
# genetic code. This writes DIR/miniprot, which calls the real miniprot with
# -T TABLE, for --miniprot_execute_path (compleasm) or the head of PATH (BUSCO).
#
# Usage: miniprot_table_wrapper.sh TABLE DIR [directory to search besides PATH ...]
set -euo pipefail
table=$1; dir=$2; shift 2
real=$(command -v miniprot || true)
for d in "$@"; do
    [ -n "$real" ] && break
    [ -x "$d/miniprot" ] && real="$d/miniprot"
done
if [ -z "$real" ]; then
    echo "miniprot_table_wrapper.sh: no miniprot in PATH${*:+ or in $*}" >&2
    exit 1
fi
# without arguments miniprot prints its options and exits with 1
usage=$("$real" 2>&1 || true)
if ! grep -q -- '-T INT' <<< "$usage"; then
    echo "miniprot_table_wrapper.sh: $real has no option -T (translation table); translation_table ${table} needs miniprot 0.13 or later" >&2
    exit 1
fi
mkdir -p "$dir"
printf '#!/bin/sh\nexec "%s" -T %s "$@"\n' "$real" "$table" > "$dir/miniprot"
chmod +x "$dir/miniprot"
