# Image of the Drusilla flow of Paludamentum (DRUSILLA_ANNOTATE, FIX_ORFS,
# GENEFINDER_LGB_FILTER): Drusilla at a pinned commit plus the packages of the
# scripts in bin/ that run next to it. The Dockerfile of the Drusilla
# repository clones Drusilla HEAD at build time and pins no package, so the
# image it gives cannot be rebuilt; this one can.
#
# The LightGBM filter (apply_lgb_model_gtf.py) loads the released text model
# with lightgbm alone. The model was converted and checked with lightgbm
# 4.7.0, which is pinned here. No scikit-learn, no joblib.
#
# Build (no build context is needed) and tag with the Drusilla version:
#   docker build -t gaiusaugustus/drusilla:0.1.0 - < docker/drusilla/Dockerfile
#   docker push gaiusaugustus/drusilla:0.1.0

# Same base as the Tiberius image: the NGC TensorFlow build is the only one
# that runs on Blackwell GPUs (sm_120).
FROM nvcr.io/nvidia/tensorflow:25.02-tf2-py3

ARG DRUSILLA_COMMIT=e3c5cf5943c05cd7c82eaccee74218804138e1a0
ARG LIGHTGBM_VERSION=4.7.0

LABEL org.opencontainers.image.title="drusilla" \
      org.opencontainers.image.description="Drusilla ORF annotator with lightgbm for the Drusilla flow of Paludamentum" \
      org.opencontainers.image.source="https://github.com/Gaius-Augustus/Drusilla" \
      drusilla.commit="${DRUSILLA_COMMIT}" \
      lightgbm.version="${LIGHTGBM_VERSION}"

USER root

# Record the Python packages of the NGC base image. The check at the end of
# this file fails the build if a pip install replaced the NGC TensorFlow or
# put PyPI CUDA/cuDNN wheels next to it (which breaks Blackwell GPUs).
RUN python3 -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name'].lower().replace('_', '-') + '==' + d.version for d in m.distributions())))" > /opt/ngc-base-packages.txt

# Drusilla's HMM layer is written against Keras 3 (tf.keras.Layer)
ENV TF_USE_LEGACY_KERAS=0
RUN python3 -m pip install --no-cache-dir --upgrade "keras>=3,<4"

# External tools Drusilla shells out to: gffread always, stringtie for BAM input
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
        git \
        gffread \
        stringtie && \
    rm -rf /var/lib/apt/lists/*

# bricks2marble[tf] and hidten[tensorflow] require tensorflow[and-cuda], whose
# and-cuda extra would install PyPI CUDA wheels on top of the NGC stack. So the
# TF-dependent packages are installed without dependencies and the rest
# explicitly (see the Tiberius Dockerfile).
RUN python3 -m pip install --no-cache-dir --no-deps bricks2marble hidten && \
    python3 -m pip install --no-cache-dir \
        numpy \
        pydantic \
        pyyaml \
        pyfaidx \
        biopython \
        pandas \
        "lightgbm==${LIGHTGBM_VERSION}"

# Drusilla at the pinned commit
RUN git init -q /opt/Drusilla && \
    cd /opt/Drusilla && \
    git fetch -q --depth 1 https://github.com/Gaius-Augustus/Drusilla ${DRUSILLA_COMMIT} && \
    git checkout -q FETCH_HEAD && \
    python3 -m pip install --no-cache-dir --no-deps .

# Model cache of the image; Paludamentum sets DRUSILLA_CACHE_DIR per task
# (params.drusilla.cache_dir) and never writes here.
ENV DRUSILLA_CACHE_DIR=/data/drusilla_cache
RUN mkdir -p /data/drusilla_cache && chmod 777 /data/drusilla_cache

# Fail the build if anything the flow needs is missing or if the NGC
# TensorFlow/CUDA stack was modified (see top of file)
RUN python3 -c "import importlib.metadata as m, sys; \
base = set(open('/opt/ngc-base-packages.txt').read().split()); \
now = {d.metadata['Name'].lower().replace('_', '-') + '==' + d.version for d in m.distributions()}; \
added = sorted(p for p in now - base if p.startswith(('nvidia-', 'tensorflow'))); \
tf = m.version('tensorflow'); \
sys.exit(f'NGC TensorFlow stack modified: tensorflow=={tf}, added/changed: {added}' if '+nv' not in tf or added else 0)" && \
    python3 -c "import lightgbm, sys; sys.exit(0 if lightgbm.__version__ == '${LIGHTGBM_VERSION}' else f'lightgbm {lightgbm.__version__} != ${LIGHTGBM_VERSION}')" && \
    python3 -c "import keras, sys; sys.exit(0 if keras.__version__.startswith('3.') else f'Keras {keras.__version__} is not Keras 3')" && \
    python3 -c "import pyfaidx, pandas, yaml, Bio" && \
    python3 -c "from drusilla.hmm.annotation_hmm import OrfAnnotationHMM" && \
    drusilla annotate --help | grep -q -- '--partial5-out' && \
    drusilla annotate --help | grep -q -- '--lorf-class' && \
    drusilla filter-subseq --help > /dev/null && \
    gffread --version && \
    stringtie --version && \
    python3 -c "import importlib.metadata as m, tensorflow as tf; b = tf.sysconfig.get_build_info(); print('drusilla', m.version('drusilla'), 'lightgbm', m.version('lightgbm'), 'keras', m.version('keras'), 'bricks2marble', m.version('bricks2marble'), 'hidten', m.version('hidten'), 'TensorFlow', tf.__version__, 'CUDA', b['cuda_version'], 'cuDNN', b['cudnn_version'])"
