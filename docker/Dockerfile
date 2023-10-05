FROM fedora:37

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gdb \
        git \
        make \
    && dnf clean all

# copy the MICM code
COPY . /micm/

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && cmake \
        -D CMAKE_BUILD_TYPE=debug \
        -D ENABLE_CLANG_TIDY:BOOL=FALSE \
        ../micm \
      && make install -j 8

WORKDIR /build