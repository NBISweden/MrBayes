# Test MrBayes+BEAGLE3 using Docker containers

- Last modified: Tue Jan 22, 2019  01:56PM
- Sign: JN

## Description

Test [MrBayes](https://github.com/NBISweden/MrBayes) and
[BEAGLE](https://github.com/beagle-dev/beagle-lib) inside
[Docker](https://www.docker.com/) containers running Linux.

The script compiles beagle and MrBayes (both serial and MPI version),
and then starts MrBayes, and tries to show the available beagle
resources (MrBayes command `showbeagle`).
Evaluation of "success" or "failure" have to be done manually
(by looking at the output).

- Note 1: the `--dns` argument is needed locally behind NRM firewall.
- Note 2: NVIDIA CUDA is not installed in the examples below.


## Test status

| Docker image | CUDA       | OpenCL | OpenMPI | Date        | Comments |
|--------------|------------|--------|-------- |-------------|----------|
| debian:sid   | Not tested | OK     | OK      | 17 Jan 2019 | Need to use `./configure LDFLAGS=-Wl,-rpath=/usr/local/lib` for beagle, or set `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"` before running `mb`|
| ubuntu:18.04 | Not tested | OK     | OK      | 17 Jan 2019 | Need to use `./configure LDFLAGS=-Wl,-rpath=/usr/local/lib` for beagle, or set `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib"` before running `mb`|
| centos:7     | Not tested | Fail   | OK      | 17 Jan 2019 | Could not get beagle to run if OpenCL was installed. Not even when using `--without-opencl` for `configure`. `showb` in MrBayes gives `OpenCL error: Unknown error from file <GPUInterfaceOpenCL.cpp>, line 115`|


## Links

- [https://en.wikipedia.org/wiki/Rpath](https://en.wikipedia.org/wiki/Rpath)
- [https://amir.rachum.com/blog/2016/09/17/shared-libraries/](https://amir.rachum.com/blog/2016/09/17/shared-libraries/)
- [https://wiki.debian.org/RpathIssue](https://wiki.debian.org/RpathIssue)
- [https://unix.stackexchange.com/questions/22926/where-do-executables-look-for-shared-objects-at-runtime](https://unix.stackexchange.com/questions/22926/where-do-executables-look-for-shared-objects-at-runtime)
- [http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html](http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html)

## Docker commands

- [centos.7.md](centos.7.md)
- [debian.sid.md](debian.sid.md)
- [ubuntu-18.04.md](ubuntu-18.04.md)
