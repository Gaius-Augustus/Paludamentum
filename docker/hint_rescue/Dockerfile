# Image of the Paludamentum hint rescue (HINT_RESCUE_LOCI, HINT_RESCUE_TIBERIUS).
#
# The hint rescue needs the Tiberius option --hints. It is on the Tiberius
# branch hint_integration, which needs bricks2marble from its branch
# intron_hints. Neither branch is part of a release. This image is the Tiberius
# image with both branches in place of the released versions, at the commits of
# the Paludamentum benchmark.
#
# The image is for the hint rescue only: its Tiberius and bricks2marble are
# older than those of the Tiberius image. It is obsolete as soon as both
# branches are released.
#
# Keep the base at larsgabriel23/tiberius:2.0.8, the image of the benchmark:
# HINT_RESCUE_LOCI needs its samtools, and this file needs its git. Both are
# gone from the Tiberius image from 2.0.9 on.
#
# Build (no build context is needed):
#   docker build -t <repository>/paludamentum-hint-rescue:<tag> - < docker/hint_rescue/Dockerfile

ARG TIBERIUS_IMAGE=larsgabriel23/tiberius:2.0.8
FROM ${TIBERIUS_IMAGE}

# Tiberius branch hint_integration and bricks2marble branch intron_hints
ARG TIBERIUS_COMMIT=5a661549e16a009d54e6aaed4a51a191177dec0f
ARG BRICKS2MARBLE_COMMIT=619c293c48c8f086cbe9fa8e1cb8b15f3f888341

LABEL org.opencontainers.image.title="paludamentum-hint-rescue" \
      org.opencontainers.image.description="Tiberius image with Tiberius hint_integration and bricks2marble intron_hints for the Paludamentum hint rescue" \
      tiberius.commit="${TIBERIUS_COMMIT}" \
      bricks2marble.commit="${BRICKS2MARBLE_COMMIT}"

USER root

# Record the Python packages of the Tiberius image. The check at the end of
# this file fails the build if anything but Tiberius and bricks2marble changed,
# so the NGC TensorFlow and CUDA stack stays as it is.
RUN python3 -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name'].lower().replace('_', '-') + '==' + d.version for d in m.distributions())))" > /opt/tiberius-image-packages.txt

# Tiberius of branch hint_integration replaces the Tiberius of the image
RUN rm -rf /opt/Tiberius && \
    git init -q /opt/Tiberius && \
    cd /opt/Tiberius && \
    git fetch -q --depth 1 https://github.com/Gaius-Augustus/Tiberius ${TIBERIUS_COMMIT} && \
    git checkout -q FETCH_HEAD && \
    python3 -m pip install --no-deps --force-reinstall . && \
    chmod +x tiberius.py tiberius/*.py tiberius/scripts/* && \
    mkdir -p model_weights && chmod -R 777 model_weights

# bricks2marble of branch intron_hints replaces the released bricks2marble
RUN python3 -m pip install --no-deps --force-reinstall \
    "bricks2marble @ git+https://github.com/Gaius-Augustus/bricks2marble@${BRICKS2MARBLE_COMMIT}"

# Fail the build if Tiberius has no --hints, if bricks2marble has no hint code,
# or if another package changed (see above)
RUN python3 /opt/Tiberius/tiberius.py --help 2>&1 | grep -q -- '--hints' && \
    python3 -c "from bricks2marble.tf.hmm.tools import get_introns_emission_distribution, get_codon_emission_distribution" && \
    python3 -c "import importlib.metadata as m, sys; \
base = set(open('/opt/tiberius-image-packages.txt').read().split()); \
now = {d.metadata['Name'].lower().replace('_', '-') + '==' + d.version for d in m.distributions()}; \
changed = sorted(p for p in now ^ base if not p.startswith(('tiberius==', 'bricks2marble=='))); \
sys.exit(f'Packages of the Tiberius image changed: {changed}' if changed else 0)" && \
    python3 -c "import importlib.metadata as m; print('tiberius', m.version('tiberius'), 'bricks2marble', m.version('bricks2marble'), 'hidten', m.version('hidten'), 'tensorflow', m.version('tensorflow'))"
