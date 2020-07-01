# Install MrBayes+Beagle3 inside an Ubuntu 20.04 Docker container

- Status (June 2020):
    - CUDA: N/A
    - OpenCL: OK
    - OpenMPI: OK
    - Beagle3: OK
    - MrBayes+MPI: OK
    - MrBayes: OK

- This file describes a "proof-of-concept" installation of the develop branch
  of MrBayes with Beagle. The installation is made inside a Docker container,
  hence, this is *not* a Dockerfile.

- Instructions for Debian are identical, but you would use a Debian container
  such as `debian:sid` instead of `ubuntu:20.04`.

- We could not get beagle to work from inside MrBayes unless Beagle is
  configured using `-rpath`, or `LD_LIBRARY_PATH` is set before running
  MrBayes.

- When configuring the beagle library, we expect to see these warning messages
  (since we asked for not using CUDA or Java):
    - `WARNING: NVIDIA CUDA nvcc compiler not found`
    - `WARNING: JDK installation not found`

- We are not building the LaTeX version of the documentation for MrBayes. For
  building the documentation, we recommend to use the `latexmk` script.  Using
  `apt install latexmk texlive-latex-extra` will install the script an all
  necessary TeX-libraries.

- When running the MPI (parallel) version of MrBayes inside the Docker
  container, we need to use `--allow-run-as-root` (since we are starting
  `mpirun` as user `root`).

## Install and run

    $ docker run -it ubuntu:20.04 /bin/bash

    # Base system
    apt update -y && apt upgrade -y

    DEBIAN_FRONTEND=noninteractive apt install -y tzdata

    apt install -y \
        autoconf \
        g++ \
        git \
        libreadline-dev \
        libtool \
        make

    # OpenCL
    apt install -y \
        ocl-icd-opencl-dev \
        pocl-opencl-icd

    # OpenMPI
    apt install -y \
        libopenmpi-dev

    # Beagle
    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    ./autogen.sh
    LDFLAGS=-Wl,-rpath=/usr/local/lib ./configure --without-jdk --disable-doxygen-doc
    make -j2
    make install

    cd /

    # MrBayes
    git clone --depth=1 --branch=develop https://github.com/NBISweden/MrBayes.git
    cd MrBayes
    ./configure --with-mpi --enable-doc=no
    make -j2

    # Test MPI (parallel) version
    mpirun --allow-run-as-root -np 1 src/mb <<< 'version;showb;quit'

    make clean
    ./configure --enable-doc=no
    make -j2

    # Test serial version
    src/mb <<< 'version;showb;quit'

