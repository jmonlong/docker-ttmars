FROM ubuntu:20.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    wget \
    gcc \ 
    samtools \
    tzdata \
    build-essential \
    bzip2 \
    git \
    sudo \
    less \
    pkg-config \
    apt-transport-https software-properties-common dirmngr gpg-agent \ 
    && rm -rf /var/lib/apt/lists/*

ENV TZ=America/Los_Angeles

WORKDIR /build

## install conda
RUN wget --quiet --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda\
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ENV PATH=/opt/conda/bin:$PATH

RUN conda update -n base -c defaults conda

RUN conda init

## LRA
RUN conda install -c bioconda lra

## samLiftOver
RUN git clone https://github.com/mchaisso/mcutils.git && \
    cd mcutils/src && make && make install
ENV PATH=/build/mcutils/bin:$PATH

## TT-Mars
RUN conda install -c anaconda numpy && \
    conda install -c conda-forge biopython

RUN conda install -c bioconda pysam && \
    conda install -c bioconda mappy && \
    conda install -c bioconda pybedtools

RUN git clone https://github.com/jmonlong/TT-Mars.git && \
    cd TT-Mars && \
    git checkout 9e8f5399243e400a9d329ca2d5ac7369c54d2a1c

WORKDIR /home
