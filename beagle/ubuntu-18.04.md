# Run MrBayes+BEAGLE3 inside Ubuntu 18.04 Docker container

- Last modified: Tue Jan 29, 2019  02:49PM
- Sign: JN
- Status: CUDA:NA, OpenCL:OK, OpenMPI:OK, BEAGLE3:OK, MB-MPI:OK, MB:OK
- Comment: Can't get beagle to work from inside mb unless beagle is configured using `-rpath`, or, `LD_LIBRARY_PATH` is set before running mb.

# Install and run

    $ docker run -it ubuntu:18.04 /bin/bash

    # Base system
    mywd=$(pwd)
    apt-get update -y
    apt-get install -y \
        libtool \
        autoconf \
        make \
        g++ \
        git \
        default-jdk \
        texlive

    # OpenCL
    apt-get install -y \
        ocl-icd-opencl-dev \
        pocl-opencl-icd

    # OpenMPI
    apt-get install -y \
        openmpi-bin \
        openmpi-common \
        libopenmpi-dev

    # BEAGLE
    cd "$mywd"
    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    ./autogen.sh
    ./configure LDFLAGS=-Wl,-rpath=/usr/local/lib
    make
    make install

    # MrBayes
    cd "$mywd"
    git clone --depth=1 --branch=develop https://github.com/NBISweden/MrBayes.git
    cd MrBayes
    ./configure --with-mpi
    make
    mpirun --allow-run-as-root -np 1 src/mb <<MBCMD
    showb
    quit
    MBCMD
    make clean
    ./configure
    make
    src/mb <<MBCMD
    showb
    quit
    MBCMD

