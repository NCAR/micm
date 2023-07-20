

Getting Started
===============

## Build and run

### CPU

### GPU

#### Checking that it compiles on your local machine

1. Build the image
```
docker build -t micm -f Dockerfile.nvhpc .
```
2. Run the container
```
docker run --rm -it micm
```
3. Compile micm. After running the previous command, you can run `make` and see your compile errors.
```
make
```
4. If you'd like, you can ssh into a running docker container and edit the files there.


## NCAR Hardware

#### On Gust
```
qinteractive -A NTDD0005 --ngpus=1
module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
mkdir build && cd build
cmake -DENABLE_OPENACC=OFF -DENABLE_CUDA=ON -D GPU_TYPE="a100" ..
make
make test
```