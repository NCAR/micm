FROM fedora:37

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc-c++ \
        gcc-fortran \
        gdb \
        git \
        make \
        netcdf-fortran-devel \
    && dnf clean all

#####
# Only needed to build the old version of micm
#####

RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

#####
#####

# copy the MICM code
COPY . /micm/

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0" \
      && cmake \
        -D CMAKE_BUILD_TYPE=debug \
        -D ENABLE_CLANG_TIDY:BOOL=FALSE \
        ../micm \
      && make 

WORKDIR /build