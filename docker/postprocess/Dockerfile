# Image of the Paludamentum post-processing: every process with the label
# 'postprocess' (compleasm, GffCompare, tRNAscan-SE, Infernal, pybarrnap,
# gene set statistics, report.html, the BUSCO lineage and Rfam downloads).
# The scripts in bin/ are not part of the image; Nextflow mounts bin/.
#
# Every tool is pinned: one micromamba environment with Bioconda and
# conda-forge packages at fixed version and build (resolved together on
# 2026-10-07), and compleasm 0.2.9 from its release archive (it bundles its
# miniprot and hmmsearch), checked against its SHA-256. pyhmmer is pinned
# as well: pybarrnap 0.5.1 asks only for pyhmmer >= 0.11 and fails on the
# keyword-only TextSequence of pyhmmer 0.12 (0.1.0 of this image: every
# pybarrnap run stopped with '__init__() takes exactly 0 positional
# arguments'); 0.5.2 needs pyhmmer >= 0.12.1. The Bioconda
# compleasm pulls Python 3.7 and HMMER 3.1b2 (through sepp), and compleasm
# 0.2.7 fails on the current file_versions.tsv of busco-data.ezlab.org.
#
# Build (no build context is needed) and tag with the version below:
#   docker build -t gaiusaugustus/paludamentum-postprocess:0.1.1 - < docker/postprocess/Dockerfile
#   docker push gaiusaugustus/paludamentum-postprocess:0.1.1

FROM mambaorg/micromamba:2.0.5

ARG COMPLEASM_VERSION=0.2.9
ARG COMPLEASM_SHA256=b57afd92ffd26b127cd2faeba5ff4eb1026e82c2f3fed0443855e0d7a5280620

LABEL org.opencontainers.image.title="paludamentum-postprocess" \
      org.opencontainers.image.description="Post-processing tools of the Paludamentum pipeline" \
      org.opencontainers.image.source="https://github.com/Gaius-Augustus/Paludamentum" \
      compleasm.version="${COMPLEASM_VERSION}"

RUN micromamba install -y -n base --override-channels -c conda-forge -c bioconda \
        python=3.12.15=h5f976f7_1_cpython \
        gffcompare=0.12.6=h9948957_4 \
        gffread=0.12.7=h077b44d_6 \
        infernal=1.1.5=pl5321h7b50bb2_4 \
        trnascan-se=2.0.12=pl5321h7b50bb2_2 \
        pybarrnap=0.5.2=pyhdfd78af_0 \
        pyhmmer=0.12.3=py312hfabe715_0 \
        matplotlib-base=3.11.2=py312hfe9aac7_1 \
        pandas=3.0.6=np2py312h91ec553_0 \
        curl=8.22.0=ha042cf0_0 \
        bzip2=1.0.8=hda65f42_10 \
        procps-ng=4.0.6=h18c060e_0 && \
    micromamba clean --all --yes

# compleasm (Python, uses pandas of the environment) with its miniprot and hmmsearch
USER root
RUN cd /tmp && \
    /opt/conda/bin/curl -fsSLO https://github.com/huangnengCSU/compleasm/releases/download/v${COMPLEASM_VERSION}/compleasm-${COMPLEASM_VERSION}_x64-linux.tar.bz2 && \
    echo "${COMPLEASM_SHA256}  compleasm-${COMPLEASM_VERSION}_x64-linux.tar.bz2" | sha256sum -c - && \
    /opt/conda/bin/bzip2 -dc compleasm-${COMPLEASM_VERSION}_x64-linux.tar.bz2 | tar --no-same-owner -x -C /opt && \
    rm compleasm-${COMPLEASM_VERSION}_x64-linux.tar.bz2 && \
    chmod -R a+rX /opt/compleasm_kit
USER $MAMBA_USER

# Nextflow runs the task scripts with bash and without the micromamba
# entrypoint, so the environment goes on PATH here
ENV PATH=/opt/conda/bin:/opt/compleasm_kit:${PATH}

# Fail the build if a tool is missing or not the pinned version
RUN set -e; \
    compleasm.py --version 2>&1 | grep -q "${COMPLEASM_VERSION}"; \
    /opt/compleasm_kit/miniprot --version > /dev/null; \
    /opt/compleasm_kit/hmmsearch -h | grep -q 'HMMER'; \
    gffcompare --version 2>&1 | grep -q '0.12.6'; \
    gffread --version | grep -qx '0.12.7'; \
    cmscan -h | grep -q 'INFERNAL 1.1.5'; \
    cmpress -h > /dev/null; \
    tRNAscan-SE -h 2>&1 | grep -q 'tRNAscan-SE 2.0.12'; \
    command -v EukHighConfidenceFilter > /dev/null; \
    pybarrnap --version | grep -q '0.5.2'; \
    pybarrnap --kingdom euk /opt/conda/lib/python3.12/site-packages/pybarrnap/example_data/fungus.fna 2> /dev/null | grep -q 'Name=28S_rRNA'; \
    curl --version > /dev/null; \
    ps --version > /dev/null; \
    python -c "import matplotlib, pandas; assert (matplotlib.__version__, pandas.__version__) == ('3.11.2', '3.0.6'), (matplotlib.__version__, pandas.__version__)"
