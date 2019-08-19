# Run MrBayes+Beagle3 inside a CentOS 7 Docker container

- Status:
    - CUDA: N/A
    - OpenCL: FAIL
    - OpenMPI: OK
    - Beagle3: OK
    - MrBayes+MPI: OK
    - MrBayes: OK

- Comment: Can't get Beagle to work from inside MrBayes unless Beagle is
configured using `-rpath`, or `LD_LIBRARY_PATH` is set before running
MrBayes.  And, can't get Beagle to run if OpenCL is installed.

# Install and run

    $ docker run -it centos:7 /bin/bash

    # Base system
    yum -y update
    yum -y install \
        autoconf \
        gcc-c++ \
        git \
        libtool \
        make \
        readline-devel

    # Note: Currently (18 Jan 2019), Beagle seems not to work if OpenCL
    # is installed. Not even when using `./configure --without-opencl`.
    # Beagle-error message when running `mb`:
    # "Unknown error from file <GPUInterfaceOpenCL.cpp>, line 115."
    #
    # OpenCL. Apparently, need to install "EPEL" sources to get OpenCL
    #yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
    #yum -y install \
    #    ocl-icd \
    #    ocl-icd-devel
    #ln -s /usr/lib64/libOpenCL.so.1 /usr/lib/libOpenCL.so # Hack https://unix.stackexchange.com/questions/292630/how-to-install-opencl-in-centos-7-using-yum

    # OpenMPI. Apparently, mpirun need to be handled by the `module` system.
    yum -y install \
        openmpi-devel
    source /etc/profile
    module load mpi

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
