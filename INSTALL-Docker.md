# Running MrBayes+Beagle3 in Docker

## Description

Run [MrBayes](https://github.com/NBISweden/MrBayes) and
[Beagle](https://github.com/beagle-dev/beagle-lib) using
[Docker](https://www.docker.com/).

The script compiles beagle and MrBayes (both serial and MPI version),
and then starts MrBayes, and tries to show the available beagle
resources (MrBayes command `showbeagle`).  Evaluation of "success" or
"failure" have to be done manually (by looking at the output).

- Note 1: A `--dns IPNUMBER` option was needed for Docker when running
in our particular local environment (firewall-related issue).
- Note 2: Nvidia CUDA was not installed in the examples below.


## Test status

| Docker image | CUDA       | OpenCL | OpenMPI | Date        | Comments |
|--------------|------------|--------|-------- |-------------|----------|
| debian:sid   | Not tested | OK     | OK      | 17 Jan 2019 | Need to use `./configure LDFLAGS=-Wl,-rpath=/usr/local/lib` for Beagle, or set `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"` before running `mb`|
| ubuntu:18.04 | Not tested | OK     | OK      | 17 Jan 2019 | Need to use `./configure LDFLAGS=-Wl,-rpath=/usr/local/lib` for Beagle, or set `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"` before running `mb`|
| centos:7     | Not tested | Fail   | OK      | 17 Jan 2019 | Could not get beagle to run if OpenCL was installed. Not even when using `--without-opencl` for `configure`. `showb` in MrBayes gives `OpenCL error: Unknown error from file <GPUInterfaceOpenCL.cpp>, line 115`|


## Docker commands

- [CentOS 7](INSTALL-Docker-CentOS.md)
- [Ubuntu 18.04 and Debian](INSTALL-Docker-Ubuntu.md)

## References

- [https://en.wikipedia.org/wiki/Rpath](https://en.wikipedia.org/wiki/Rpath)
- [https://amir.rachum.com/blog/2016/09/17/shared-libraries/](https://amir.rachum.com/blog/2016/09/17/shared-libraries/)
- [https://wiki.debian.org/RpathIssue](https://wiki.debian.org/RpathIssue)
- [https://unix.stackexchange.com/questions/22926/where-do-executables-look-for-shared-objects-at-runtime](https://unix.stackexchange.com/questions/22926/where-do-executables-look-for-shared-objects-at-runtime)
- [http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html](http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html)
