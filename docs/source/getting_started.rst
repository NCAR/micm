

Getting Started
===============

Build and Run
-------------

CPU
~~~
To build and install MICM locally, you must have the following libraries installed:

- [sphinx](https://github.com/sphinx-doc/sphinx)
- [sphinx-book-theme](https://github.com/executablebooks/sphinx-book-theme)
- [sphinx-design](https://github.com/executablebooks/sphinx-design)
- [breathe](https://github.com/breathe-doc/breathe)

You must also have CMake installed on your machine.

Open a terminal window, navigate to a folder where you would like the MICM files to exist,
and run the following commands:

```
git clone https://github.com/NCAR/micm.git
cd micm
mkdir build
cd build
ccmake ..
sudo make install -j 8
```

GPU
~~~

Checking that it compiles on your local machine

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


NCAR Hardware
-------------

On Gust
~~~~~~~
```
qinteractive -A NTDD0005 --ngpus=1
module load cmake/3.25.2 nvhpc/23.1 cuda/11.7.1
mkdir build && cd build
cmake -DENABLE_OPENACC=OFF -DENABLE_CUDA=ON -D GPU_TYPE="a100" ..
make
make test
```
