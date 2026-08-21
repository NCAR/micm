FROM fedora:latest

ARG BUILD_TYPE=release

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gdb \
        ninja-build \
        git \
        make \
    && dnf clean all

# copy the MICM code
COPY . /micm/

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && cmake \
        -G Ninja \
        -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -D MICM_ENABLE_CLANG_TIDY:BOOL=FALSE \
        ../micm \
      && ninja install -j 8

WORKDIR /build