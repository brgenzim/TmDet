# Use Ubuntu 24.04 as the base image
FROM ubuntu:24.04

ARG GEMMI_VERSION=0.7.0

# Update the package list and install essential build tools and libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    git \
    libcurl4-openssl-dev \
    libzip-dev \
    nlohmann-json3-dev

# Install the latest GCC compiler that supports the C++20 standard
RUN apt-get install -y gcc g++ gdb

# Install Eigen3 from the Ubuntu repositories
RUN apt-get install -y libeigen3-dev

# Download, build, and install Gemmi
WORKDIR /opt
RUN curl -L -O https://github.com/project-gemmi/gemmi/archive/refs/tags/v$GEMMI_VERSION.tar.gz && \
    tar -xzf v$GEMMI_VERSION.tar.gz && \
    rm v$GEMMI_VERSION.tar.gz && \
    cd gemmi-$GEMMI_VERSION && \
    mkdir build && \
    cmake -B build && \
    make -C build && \
    make -C build install

# Download pugixml 1.14
WORKDIR /usr/local/src/contrib
RUN curl -L -O https://github.com/zeux/pugixml/archive/refs/tags/v1.14.tar.gz && \
    tar -xzf v1.14.tar.gz && \
    rm v1.14.tar.gz && \
    ln -s pugixml-1.14 pugixml

RUN addgroup --gid 1120 tusnady

# Set the working directory for tmdet source and build
WORKDIR /usr/local/src/tmdet

# Default command (can be used to compile and run projects)
CMD ["/bin/bash"]
