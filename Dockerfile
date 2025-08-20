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

RUN meson install -C project/builddir --destdir /AppDir

RUN wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage -O /linuxdeploy.AppImage && \
    chmod +x /linuxdeploy.AppImage && \
    /linuxdeploy.AppImage --appimage-extract && \
    OUTPUT="Precision_Mesh.AppImage" ./squashfs-root/AppRun \
        --appdir AppDir \
        --desktop-file /AppDir/usr/share/applications/precision_mesh.desktop \
        --output appimage && \
    rm /linuxdeploy.AppImage && \
    rm -rf ./squashfs-root


FROM ubuntu:22.04
COPY --from=builder /project/builddir/precision_mesh /usr/bin/precision_mesh

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libfontconfig1 \
        libmpfr6 \
        libtbb12 \
        libx11-6 \
    && rm -rf /var/lib/apt/lists/*

CMD precision_mesh