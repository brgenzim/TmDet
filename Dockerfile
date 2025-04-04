# Use Ubuntu 24.04 as the base image
FROM ubuntu:24.04

ARG GEMMI_VERSION=0.7.0
ARG PUGIXML_VERSION=1.14
ARG SRC_DESTINATION=/usr/local/src/tmdet
ARG ENV_FILE=/etc/tmdet.env

# Update the package list and install essential build tools and libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    git \
    libcurl4-openssl-dev \
    libeigen3-dev \
    libzip-dev \
    nlohmann-json3-dev

# Install the latest GCC compiler that supports the C++20 standard
RUN apt-get install -y gcc g++ gdb

# Download, build, and install Gemmi
WORKDIR /opt
RUN curl -L -O https://github.com/project-gemmi/gemmi/archive/refs/tags/v$GEMMI_VERSION.tar.gz && \
    tar -xzf v$GEMMI_VERSION.tar.gz && \
    rm v$GEMMI_VERSION.tar.gz && \
    cd gemmi-$GEMMI_VERSION && \
    mkdir build && \
    cmake -B build && \
    make -j4 -C build && \
    make -C build install

# Download pugixml
WORKDIR /usr/local/src/contrib
RUN curl -L -O https://github.com/zeux/pugixml/archive/refs/tags/v$PUGIXML_VERSION.tar.gz && \
    tar -xzf v$PUGIXML_VERSION.tar.gz && \
    rm v$PUGIXML_VERSION.tar.gz && \
    ln -s pugixml-$PUGIXML_VERSION pugixml

# Set the working directory for tmdet source and build
COPY . $SRC_DESTINATION
WORKDIR $SRC_DESTINATION
RUN cmake -B build && \
    make -j4 -C build && \
    make -C build install && \
    cp .env.example $ENV_FILE && \
    chmod 644 $ENV_FILE

RUN rm -rf $SRC_DESTINATION

# Default command to run tmdet
ENV LD_LIBRARY_PATH=/usr/local/lib
WORKDIR /work
ENTRYPOINT [ "/usr/local/bin/tmdet", "-e", "/etc/tmdet.env" ]
CMD [ "-h" ]
