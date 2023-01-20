FROM fedora:35

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gcc-fortran \
        git \
        make \
    && dnf clean all

# copy the MICM code
COPY . /micm/

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && cmake ../micm \
      && make -j 8

WORKDIR /build