FROM fedora:34

RUN dnf -y update \
    && dnf -y install \
        gcc-c++ \
        gcc \
        gdb \
        cmake \
        make \
        lcov \
        valgrind \
    && dnf clean all

COPY . /micm/
RUN mkdir /build \
      && cd /build \
      && cmake /micm \
      && make
