FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        file \
        git \
        libboost-dev \
        libfontconfig-dev \
        libgmp3-dev \
        libmpfr-dev \
        libtbb-dev \
        libx11-dev \
        meson \
        wget \
    && rm -rf /var/lib/apt/lists/*

COPY . project

RUN cd project && \
    meson setup --prefix=/usr builddir

RUN cd project/builddir && \
    meson compile

FROM ubuntu:22.04
COPY --from=builder /project/builddir/precision_mesh /usr/bin/precision_mesh

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libfontconfig1 \
        libmpfr6 \
        libtbb12 \
        libx11-6 \
    && rm -rf /var/lib/apt/lists/*

CMD ["precision_mesh"]