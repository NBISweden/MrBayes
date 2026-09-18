# Run MrBayes+Beagle3 inside a CentOS 8 Docker container

- Status (June 2020):
    - CUDA: N/A
    - OpenCL: FAIL
    - OpenMPI: OK
    - Beagle3: OK
    - MrBayes+MPI: OK
    - MrBayes: OK

- This file describes a "proof-of-concept" installation of the develop branch
  of MrBayes with Beagle. The installation is made inside a Docker container,
  hence, this is *not* a Dockerfile.

- We could not get beagle to work from inside MrBayes unless Beagle is
  configured using `-rpath`, or `LD_LIBRARY_PATH` is set before running
  MrBayes.

- When configuring the beagle library, we expect to see these warning messages
  (since we asked for not using OpenCL, CUDA or Java):
    - `WARNING: OpenCL not found or disabled`
    - `WARNING: NVIDIA CUDA nvcc compiler not found`
    - `WARNING: JDK installation not found`

- We are not building the LaTeX version of the documentation for MrBayes. For
  building the documentation, we recommend to use the `latexmk` script. Using
  the commands below will install the script an all necessary TeX-libraries.

        yum install -y epel-release
        yum -y update
        yum -y install texlive latexmk texlive-wrapfig

- When running the MPI (parallel) version of MrBayes inside the Docker
  container, we need to use `--allow-run-as-root` (since we are starting
  `mpirun` as user `root`).


# Install and run

    $ docker run -it centos:8 /bin/bash

    LC_ALL=C

    # Base system
    yum -y install \
        autoconf \
        gcc-c++ \
        git \
        libtool \
        make \
        readline-devel

    # OpenCL
    # Tue 30 Jun 2020: Beagle compiles with OpenCL, but when run in MrBayes, we get
    # OpenCL error: Unknown error from file <GPUInterfaceOpenCL.cpp>, line 115.
    #yum -y install \
    #    'dnf-command(config-manager)'
    #yum -y config-manager --set-enabled PowerTools
    #yum -y update
    #yum -y install \
    #    ocl-icd \
    #    ocl-icd-devel
    #ln -s /usr/lib64/libOpenCL.so.1 /usr/lib/libOpenCL.so # https://unix.stackexchange.com/questions/292630/how-to-install-opencl-in-centos-7-using-yum

    # OpenMPI
    yum -y install \
        openmpi-devel
    source /etc/profile
    module load mpi

    # Beagle
    git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
    cd beagle-lib
    ./autogen.sh
    LDFLAGS=-Wl,-rpath=/usr/local/lib ./configure --without-jdk --without-cuda --without-opencl --disable-doxygen-doc
    make -j2
    make install

    cd /

    # MrBayes
    git clone --depth=1 --branch=develop https://github.com/NBISweden/MrBayes.git
    cd MrBayes
    ./configure --with-mpi --enable-doc=no
    make -j2

    mpirun --allow-run-as-root -np 1 src/mb <<< 'version;showb;quit'

    make clean
    ./configure --enable-doc=no
    make -j2

    src/mb <<< 'version;showb;quit'

