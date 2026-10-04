# Image of the Paludamentum evidence tools: every process with the label
# 'container' (StringTie, HISAT2, minimap2, miniprot, miniprothint, DIAMOND,
# TransDecoder, TD2, bam2hints, SRA downloads, the scripts in bin/). Tiberius runs
# in the Tiberius image (label 'tiberius'), Vipsania, Drusilla and the hint
# rescue in theirs.
#
# Every tool is pinned. The versions are those of the image the pipeline was
# benchmarked with (larsgabriel23/tiberius:2.0.8, which cloned Augustus,
# miniprot, miniprothint and miniprot-boundary-scorer at their HEAD), so the
# outputs do not change. Ubuntu 24.04 is the base of that image as well; the
# apt packages are pinned to its versions.
#
# Build (no build context is needed) and tag with the version below:
#   docker build -t gaiusaugustus/paludamentum-evidence:0.2.1 - < docker/evidence/Dockerfile
#   docker push gaiusaugustus/paludamentum-evidence:0.2.1
#
# 0.2.1: TD2 1.1.0 (the ORF finder of the TransDecoder flow by default, see
# conf/hc_genes.yaml) with psauron 1.1.3 and the CPU build of PyTorch 2.14.0,
# the versions of docs/orf_finder_comparison.md. (0.2.0 lacked psauron on
# PATH, so TD2.Predict failed; do not use it.)

# ---------------------------------------------------------------- build stage
FROM ubuntu:24.04 AS build

ARG AUGUSTUS_COMMIT=220e5a63c0fd546563472ec6d1f9271faddaf569
ARG MINIPROT_COMMIT=81f9b93481cec16622eff339fc998e4f38da344e
ARG MINIPROT_BOUNDARY_SCORER_COMMIT=b0d103ef920dc11a567700c1b3dec022fea8de6d
ARG MINIPROTHINT_COMMIT=07cc5abe7fb83d4cbf98a3342e8423beefa81c99

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        bzip2 \
        zlib1g-dev \
        libbamtools-dev=2.5.2+dfsg-4 && \
    rm -rf /var/lib/apt/lists/*

# fetch <repository> <commit> <directory>
RUN printf '%s\n' '#!/bin/sh' 'set -e' \
        'git init -q "$3"' \
        'git -C "$3" fetch -q --depth 1 "$1" "$2"' \
        'git -C "$3" checkout -q FETCH_HEAD' > /usr/local/bin/fetch && \
    chmod +x /usr/local/bin/fetch

RUN fetch https://github.com/lh3/miniprot ${MINIPROT_COMMIT} /src/miniprot && \
    make -C /src/miniprot && \
    mkdir -p /opt/miniprot && cp /src/miniprot/miniprot /opt/miniprot/

RUN fetch https://github.com/tomasbruna/miniprot-boundary-scorer ${MINIPROT_BOUNDARY_SCORER_COMMIT} /src/mbs && \
    make -C /src/mbs && \
    mkdir -p /opt/miniprot-boundary-scorer && cp /src/mbs/miniprot_boundary_scorer /opt/miniprot-boundary-scorer/

# Only bam2hints of Augustus (BAM2HINTS); nothing else of Augustus is used
RUN fetch https://github.com/Gaius-Augustus/Augustus ${AUGUSTUS_COMMIT} /src/augustus && \
    make -C /src/augustus/auxprogs/bam2hints && \
    mkdir -p /opt/augustus/bin && cp /src/augustus/auxprogs/bam2hints/bam2hints /opt/augustus/bin/

# miniprothint is Python (standard library only); it calls its own scripts
RUN fetch https://github.com/tomasbruna/miniprothint ${MINIPROTHINT_COMMIT} /opt/miniprothint && \
    rm -rf /opt/miniprothint/.git /opt/miniprothint/example

# Release archives, checked against their SHA-256
RUN mkdir -p /dl && cd /dl && \
    curl -fsSLO https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz && \
    curl -fsSLO https://github.com/gpertea/stringtie/releases/download/v3.0.3/stringtie-3.0.3.Linux_x86_64.tar.gz && \
    curl -fsSLO https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 && \
    curl -fsSLO https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    curl -fsSLO https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.3.0/sratoolkit.3.3.0-ubuntu64.tar.gz && \
    printf '%s\n' \
        '41dd5e95f6ba946ff21340417d867e5e99f123b4035779b25d3cffd20b828a30  TransDecoder-v5.7.1.tar.gz' \
        'dcd8b0c0069d0707db2decc9ecb7e7455b4e3b3f6c08ef13744860a9fb03acef  stringtie-3.0.3.Linux_x86_64.tar.gz' \
        'b4c0677f282f448dac7bc1cd2fe5a45d12cb07a8b1edc316a45a26bbf75b7f07  minimap2-2.30_x64-linux.tar.bz2' \
        '56ed696812cde1b5a08556d66fbd48cd0f2037cea56907500861bba291245005  gffread-0.12.7.Linux_x86_64.tar.gz' \
        'b84fceeeca3ccfab2eac88dc860a1a71c830e1dd95cea9cc53e2f93668e7dcc0  sratoolkit.3.3.0-ubuntu64.tar.gz' \
        | sha256sum -c - && \
    tar --no-same-owner -xzf TransDecoder-v5.7.1.tar.gz -C /opt && \
    tar --no-same-owner -xzf stringtie-3.0.3.Linux_x86_64.tar.gz -C /opt && \
    tar --no-same-owner -xjf minimap2-2.30_x64-linux.tar.bz2 -C /opt && \
    tar --no-same-owner -xzf gffread-0.12.7.Linux_x86_64.tar.gz -C /opt && \
    tar --no-same-owner -xzf sratoolkit.3.3.0-ubuntu64.tar.gz -C /opt

# ---------------------------------------------------------------- final image
FROM ubuntu:24.04

ARG AUGUSTUS_COMMIT=220e5a63c0fd546563472ec6d1f9271faddaf569
ARG MINIPROT_COMMIT=81f9b93481cec16622eff339fc998e4f38da344e
ARG MINIPROT_BOUNDARY_SCORER_COMMIT=b0d103ef920dc11a567700c1b3dec022fea8de6d
ARG MINIPROTHINT_COMMIT=07cc5abe7fb83d4cbf98a3342e8423beefa81c99

LABEL org.opencontainers.image.title="paludamentum-evidence" \
      org.opencontainers.image.description="Evidence tools of the Paludamentum pipeline" \
      org.opencontainers.image.source="https://github.com/Gaius-Augustus/Paludamentum" \
      augustus.commit="${AUGUSTUS_COMMIT}" \
      miniprot.commit="${MINIPROT_COMMIT}" \
      miniprot-boundary-scorer.commit="${MINIPROT_BOUNDARY_SCORER_COMMIT}" \
      miniprothint.commit="${MINIPROTHINT_COMMIT}"

# samtools, HISAT2 and DIAMOND from Ubuntu (the versions that the benchmark
# image ran: its DIAMOND 2.1.16 tarball came after the apt DIAMOND on PATH).
# libbamtools for bam2hints, liburi-perl for TransDecoder, procps for the
# Nextflow task metrics.
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
        samtools=1.19.2-1build2 \
        hisat2=2.2.1-4 \
        diamond-aligner=2.1.9-1 \
        libbamtools-dev=2.5.2+dfsg-4 \
        liburi-perl=5.27-1 \
        perl \
        python3 \
        python3-venv \
        procps \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Python packages of the scripts in bin/, in a venv first on PATH (Ubuntu's
# Python is externally managed); same versions as in the benchmark image
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir \
        "biopython==1.88" \
        "numpy==1.26.4" \
        "pandas==2.2.3" \
        "pyyaml==6.0.2"

# TD2 in its own venv, so that its dependencies do not change the packages of
# the bin/ scripts; only TD2.LongOrfs, TD2.Predict and psauron (which
# TD2.Predict runs as a command) go on PATH
RUN python3 -m venv /opt/td2-venv && \
    /opt/td2-venv/bin/pip install --no-cache-dir \
        --index-url https://download.pytorch.org/whl/cpu \
        --extra-index-url https://pypi.org/simple \
        "torch==2.14.0+cpu" \
        "psauron==1.1.3" \
        "TD2==1.1.0" && \
    mkdir -p /opt/td2/bin && \
    ln -s /opt/td2-venv/bin/TD2.LongOrfs /opt/td2-venv/bin/TD2.Predict /opt/td2-venv/bin/psauron /opt/td2/bin/

COPY --from=build /opt/ /opt/

ENV PATH=/opt/venv/bin:/opt/td2/bin:/opt/miniprot:/opt/miniprot-boundary-scorer:/opt/miniprothint:/opt/augustus/bin:/opt/TransDecoder-TransDecoder-v5.7.1:/opt/TransDecoder-TransDecoder-v5.7.1/util:/opt/stringtie-3.0.3.Linux_x86_64:/opt/minimap2-2.30_x64-linux:/opt/gffread-0.12.7.Linux_x86_64:/opt/sratoolkit.3.3.0-ubuntu64/bin:${PATH}

# Fail the build if a tool is missing or not the pinned version
RUN set -e; \
    samtools --version | head -1 | grep -qx 'samtools 1.19.2'; \
    hisat2 --version | grep -q 'version 2.2.1'; \
    hisat2-build --version > /dev/null; \
    diamond version | grep -qx 'diamond version 2.1.9'; \
    stringtie --version | grep -qx '3.0.3'; \
    minimap2 --version | grep -qx '2.30-r1287'; \
    gffread --version | grep -qx '0.12.7'; \
    fasterq-dump --version | grep -q '3.3.0'; \
    miniprot --version | grep -qx '0.18-r281'; \
    test -x /opt/miniprot-boundary-scorer/miniprot_boundary_scorer; \
    miniprothint.py --help > /dev/null; \
    test -x /opt/augustus/bin/bam2hints; \
    TransDecoder.LongOrfs --version | grep -q '5.7.1'; \
    /opt/td2-venv/bin/python -c "from importlib.metadata import version as v; assert (v('td2'), v('psauron'), v('torch')) == ('1.1.0', '1.1.3', '2.14.0+cpu'), (v('td2'), v('psauron'), v('torch'))"; \
    cd /tmp; \
    python3 -c "import random; random.seed(1); c = [a+b+d for a in 'ACGT' for b in 'ACGT' for d in 'ACGT' if a+b+d not in ('TAA', 'TAG', 'TGA')]; u = lambda n: ''.join(random.choice('ACGT') for _ in range(n)); o = 'ATG' + ''.join(random.choice(c) for _ in range(300)) + 'TAA'; open('td2.fa', 'w').write('>t1\\n' + u(120) + o + u(150) + '\\n')"; \
    TD2.LongOrfs -t td2.fa -O td2.dir > /dev/null; \
    TD2.Predict -t td2.fa -O td2.dir > /dev/null; \
    grep -q 'ORF type:complete (+).* t1:121-1026(+)' td2.fa.TD2.pep; \
    rm -rf /tmp/td2.*; cd /; \
    for s in gtf_genome_to_cdna_fasta.pl cdna_alignment_orf_to_genome_orf.pl gtf_to_alignment_gff3.pl; do command -v $s > /dev/null; done; \
    perl -MURI::Escape -e 1; \
    python3 -c "import Bio, numpy, pandas, yaml; print('biopython', Bio.__version__, 'numpy', numpy.__version__, 'pandas', pandas.__version__, 'pyyaml', yaml.__version__)"
