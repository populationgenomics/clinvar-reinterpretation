FROM python:3.10-bullseye

RUN apt update && apt install -y \
        apt-transport-https \
        bzip2 \
        ca-certificates \
        git \
        gnupg \
        openjdk-11-jdk-headless \
        wget \
        zip && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

# install bcftools
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git \
    && git clone https://github.com/samtools/bcftools.git \
    && cd bcftools \
    && make \
    && strip bcftools plugins/*.so \
    && make install \
    && cd - \
    && rm -r bcftools htslib

# now do some fun stuff, installing ClinvArbitration
WORKDIR /clinvarbitration

COPY example_script_docker.sh pyproject.toml README.md ./
COPY src src/
COPY bcftools_data bcftools_data/

RUN pip install .
