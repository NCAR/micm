FROM fedora:29

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        netcdf-fortran-devel \
        cmake \
        make \
    && dnf clean all

# copy the MusicBox code
COPY . /musica-core/

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.1.0.tar.gz \
    && tar -zxvf 8.1.0.tar.gz \
    && cd json-fortran-8.1.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

# get a tag and build the model
RUN mkdir /build \
    && cd /build \
    && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.1.0" \
    && cmake ../musica-core \
    && make
