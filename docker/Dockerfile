FROM fedora:37

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gdb \
        git \
        make \
        json-devel \
    && dnf clean all

# copy the MICM code
COPY . /micm/

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && cmake \
        -D CMAKE_BUILD_TYPE=release \
        -D MICM_ENABLE_CLANG_TIDY:BOOL=FALSE \
        ../micm \
      && make install -j 8

# now test if we can use the installed files
RUN cd /micm/test/integration/cmake/find_package \
      && mkdir build && cd build \
      && cmake .. \
      && make \
      && make test

# now test that fetch content can be used
RUN cd /micm/test/integration/cmake/fetch_content \
      && mkdir build && cd build \
      && cmake .. \
      && make \
      && make test

WORKDIR /build