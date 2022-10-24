FROM fedora:29

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        netcdf-fortran-devel \
        cmake \
        make \
    && dnf clean all

# copy the MICM code
COPY . /micm/

# use the test preprocessor output
RUN cp /micm/test/preprocessor_output/* /micm/src/preprocessor_output/

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

# build the library and run the tests
RUN mkdir /build \
      && cd /build \
      && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0" \
      && cmake -D ENABLE_UTIL_ONLY=ON \
               ../micm \
      && make
