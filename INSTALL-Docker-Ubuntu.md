# Run MrBayes+Beagle3 inside an Ubuntu 18.04 Docker container

- Status:
    - CUDA: N/A
    - OpenCL: OK
    - OpenMPI: OK
    - Beagle3: OK
    - MrBayes+MPI: OK
    - MrBayes: OK

- Comment: Can't get beagle to work from inside MrBayes unless Beagle is
configured using `-rpath`, or `LD_LIBRARY_PATH` is set before running
MrBayes.

- Instructions for Debian are identical, but you would use a Debian
container such as `debian:sid` instead of `ubuntu:18.04`.

# Install and run

    $ docker run -it ubuntu:18.04 /bin/bash

    # Base system
    apt-get update -y
    apt-get install -y \
        autoconf \
        g++ \
        git \
        libreadline-dev \
        libtool \
        make

    # OpenCL
    apt-get install -y \
        ocl-icd-opencl-dev \
        pocl-opencl-icd

    # OpenMPI
    apt-get install -y \
        libopenmpi-dev

    # Beagle
    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    ./autogen.sh
    ./configure --without-jdk LDFLAGS=-Wl,-rpath=/usr/local/lib
    make -j2
    make install
    cd ..

    # MrBayes
    git clone --depth=1 --branch=develop https://github.com/NBISweden/MrBayes.git
    cd MrBayes
    ./configure --with-mpi
    make -j2

    mpirun --allow-run-as-root -np 1 src/mb <<END_MRBAYES
    showb
    quit
    END_MRBAYES

    make clean
    ./configure
    make -j2

    src/mb <<END_MRBAYES
    showb
    quit
    END_MRBAYES
