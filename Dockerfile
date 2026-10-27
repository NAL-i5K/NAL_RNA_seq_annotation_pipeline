# Pin an LTS release instead of "latest" so rebuilds are reproducible
FROM ubuntu:24.04

ENV PIP_DISABLE_PIP_VERSION_CHECK=1

ARG RELEASE_VERSION=2026.10.26
LABEL org.opencontainers.image.title="NAL RNA-Seq Annotation Pipeline" \
      org.opencontainers.image.version="${RELEASE_VERSION}" \
      org.opencontainers.image.source="https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline"

ENV DEBIAN_FRONTEND=noninteractive

ARG SAMTOOLS_VERSION=1.24
ARG SRATOOLKIT_VERSION=3.4.1

# System packages in one layer; apt cache cleaned to reduce image size
RUN apt-get update && apt-get install -y \
      ca-certificates \
      rsync \
      nano vim \
      python3 python3-venv python3-dev python-is-python3 \
      perl libxml-libxml-perl liburi-perl \
      default-jre-headless \
      git cmake \
      wget bzip2 make gcc \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses-dev \
      libxml2-dev \
      ncbi-entrez-direct \
      rsem \
 && rm -rf /var/lib/apt/lists/*

# samtools: use the release tarball, which bundles HTSlib
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && cd samtools-${SAMTOOLS_VERSION} \
 && ./configure --prefix=/usr/local \
 && make -j"$(nproc)" \
 && make install \
 && cd / && rm -rf samtools-${SAMTOOLS_VERSION}*

# SRA Toolkit
RUN wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-alma_linux64.tar.gz \
 && tar -xzf sratoolkit.${SRATOOLKIT_VERSION}-alma_linux64.tar.gz -C /opt \
 && rm sratoolkit.${SRATOOLKIT_VERSION}-alma_linux64.tar.gz
ENV PATH="/opt/sratoolkit.${SRATOOLKIT_VERSION}-alma_linux64/bin:${PATH}"

# wigToBigWig (HTTPS, sensible permissions)
#RUN wget -q -O /usr/local/bin/wigToBigWig \
#      https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig \
# && chmod 755 /usr/local/bin/wigToBigWig

RUN rsync -a rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/wigToBigWig /usr/local/bin/wigToBigWig \
 && chmod 755 /usr/local/bin/wigToBigWig \
 && (/usr/local/bin/wigToBigWig 2>&1 || true) | grep -q "wigToBigWig v"

# Python virtual environment (avoids PEP 668 "externally-managed-environment")
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"
RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
 && pip install --no-cache-dir pysam

# pysam must still provide the deprecated names the script uses
RUN python3 -c "import pysam; pysam.Samfile; pysam.index; print('pysam', pysam.__version__)"

# Running the script with no args imports it (exercising the signal-table parsing)
# and runs its dependency check. Fail the build if anything is missing.
RUN out="$(python3 ${RNANNOT_LIB}/bam_to_bigwig/bam_to_bigwig.py 2>&1)"; \
    echo "$out" | tail -5; \
    ! echo "$out" | grep -q -e "Missing dependency" -e "Traceback"

# Install the pipeline
RUN mkdir -p /opt/RNA_repo
ARG FASTQC_VERSION=0.12.1
ARG TRIMMOMATIC_VERSION=0.39
ARG HISAT2_VERSION=2.2.1
ARG BBMAP_VERSION=39.33
ARG PICARD_VERSION=3.4.0
ARG REGTOOLS_VERSION=1.0.0

ENV RNANNOT_LIB=/opt/RNA_repo/rnannot/lib
WORKDIR ${RNANNOT_LIB}

# FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
 && unzip -q fastqc_v${FASTQC_VERSION}.zip && rm fastqc_v${FASTQC_VERSION}.zip \
 && chmod 755 FastQC/fastqc

# Trimmomatic (verify the URL; the usadellab site has moved releases to GitHub)
RUN wget -q https://github.com/usadellab/Trimmomatic/releases/download/v${TRIMMOMATIC_VERSION}/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
 && unzip -q Trimmomatic-${TRIMMOMATIC_VERSION}.zip && rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip

# HISAT2 2.2.1 prebuilt binary
RUN wget -q -O hisat2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
 && unzip -q hisat2.zip && rm hisat2.zip \
 && chmod 755 hisat2-${HISAT2_VERSION}/hisat2*

# BBMap
RUN wget -q -O bbmap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download \
 && tar -xzf bbmap.tar.gz && rm bbmap.tar.gz

# Picard
RUN wget -q -O picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar

# regtools (pinned tag; confirm the tag name with: git ls-remote --tags https://github.com/griffithlab/regtools)
RUN git clone --depth 1 --branch ${REGTOOLS_VERSION} https://github.com/griffithlab/regtools \
 && mkdir regtools/build && cd regtools/build && cmake .. && make -j"$(nproc)"

# bam_to_bigwig (replace <commit> with a specific commit hash)
RUN git clone https://github.com/NAL-i5K/bam_to_bigwig.git \
 && cd bam_to_bigwig && git checkout 715778648040763912e7fe048a8957d77cb6a1ff
COPY . /opt/RNA_repo/
WORKDIR /opt/RNA_repo
RUN pip install --no-cache-dir -e .

# Output folder for bind-mounting
RUN mkdir -p /opt/output
WORKDIR /opt/output
