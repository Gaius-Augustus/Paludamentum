#!/usr/bin/env bash
# The GPU that FANTASIA-Lite needs: a visible CUDA device with at least
# 15 GB free memory (ProtT5-XL uses about 14 GB). Run by FANTASIA_GPU_CHECK
# at the start of the run and again by FANTASIA_ANNOTATE before the model
# loads (modules/fantasia.nf).
set -euo pipefail

if ! command -v nvidia-smi > /dev/null 2>&1 || ! nvidia-smi -L > /dev/null 2>&1; then
    echo "FANTASIA-Lite needs a CUDA GPU, but no GPU is visible on $(hostname). Give the label gpu a GPU queue in your config, or set fantasia.run = false." >&2
    exit 1
fi
gpu=${CUDA_VISIBLE_DEVICES:-0}
gpu=${gpu%%,*}
free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i "$gpu" 2>/dev/null | tr -d ' ' || echo 0)
if [ "${free:-0}" -lt 15000 ] 2>/dev/null; then
    echo "GPU $gpu has $free MiB free; ProtT5-XL needs at least 15000 MiB." >&2
    exit 1
fi
