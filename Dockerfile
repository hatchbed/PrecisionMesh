FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libboost-dev \
        libfontconfig-dev \
        libgmp3-dev \
        libmpfr-dev \
        libtbb-dev \
        meson \
    && rm -rf /var/lib/apt/lists/*

COPY . project

RUN cd project && \
    meson setup builddir

RUN cd project/builddir && \
    meson compile
