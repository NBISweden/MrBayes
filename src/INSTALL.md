MrBayes 3.2.7-dev (development not-yet-3.2.7 version)
========================================================================

Building MrBayes using the GNU autotools (Autoconf and Automake) from a bare-bones Git clone:
------------------------------------------------------------------------

The code in the Git repository for MrBayes is (as of writing) set up
to build the software using the GNU autotools, i.e. using Automake and
Autoconf, on a Unix-like system.  Therefore,

* Make sure [Automake 1.15](http://www.gnu.org/software/automake/)
  is installed.

* Make sure [Autoconf 2.69](http://www.gnu.org/software/autoconf/)
  is installed.

The `configure.ac` script here uses the `pkg-config` tool for detection
of external libraries (Beagle) and also uses m4 macros from the
autoconf-archive collection for detection of MPI, Readline, and some
architecture-dependent compiler flags.  Therefore,

* Make sure
  [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
  is installed.

* Make sure
  [autoconf-archive 2016.03.20](http://www.gnu.org/software/autoconf-archive/)
  (or later) is installed.

* On some systems that uses metaauto for wrapping multiple versions of
  the GNU autotools (OpenBSD for example):

        $ export AUTOCONF_VERSION=2.69
        $ export AUTOMAKE_VERSION=1.15

When all the prerequisites are in place:

* Run `autoreconf` to generate the `configure` script:

        $ autoreconf -i

* Run `configure` to create `Makefile` (see `./configure --help` for
  info about enabling and disabling features of MrBayes, or how to
  install in non-standard locations etc.):

        $ ./configure

* Make the project:

        $ make

* Install it on the system (optional):

        $ make install

Picking up the Beagle library:
------------------------------------------------------------------------

Assuming [Beagle](https://github.com/beagle-dev/beagle-lib) was compiled
and installed in a standard location, it will be picked up by the
`configure` script automatically (using the `pkg-config` tool).  If it
isn't found, or the `configure` script is run with `--without-beagle`,
Beagle will not be used.

If Beagle is installed in a non-default location, point the
`PKG_CONFIG_PATH` environment variable to the `pkgconfig` directory that
Beagle installed its `hmsbeagle-1.pc` file into.

For example, I installed Beagle in my home directory using

    $ ./configure --prefix=$HOME/local
    $ make && make install

This means I will need to point `PKG_CONFIG_PATH` to
`$HOME/local/lib/pkgconfig` when I build MrBayes:

    $ env PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig ./configure
    [...]
    checking for pkg-config... /usr/bin/pkg-config
    checking pkg-config is at least version 0.9.0... yes
    checking for BEAGLE... yes
    [...]
    $ make && make install

It's also possible to set the `BEAGLE_CFLAGS` and `BEAGLE_LIBS`
environment variables prior to running `configure`:

    $ env BEAGLE_CFLAGS="-I$HOME/local/include/hmsbeagle-1" \
          BEAGLE_LIBS="-L$HOME/local/lib -lhmsbeagle" \
          ./configure

Picking up and using MPI:
------------------------------------------------------------------------

By default, the `configure` script will find any available MPI C
compiler and use it.  This behaviour may be changed by using the
`--without-mpi` flag when running the `configure`.  Setting the
environment variable `CC` to a specific C compiler will also disable MPI
support (unless this compiler happens to be an MPI C compiler).

When running `configure`, unless `--without-mpi` is used, it will always
say it's going to use MPI, even though no MPI C compiler is available:

    $ ./configure
    [...]
    checking whether to compile using MPI... yes
    [...]

Don't worry, MPI support will not actually be enabled if configure can't
find the MPI libraries and header files:

    [...]
    checking for function MPI_Init... no
    checking for function MPI_Init in -lmpi... no
    checking for function MPI_Init in -lmpich... no
    [...]

You may also set the `MPICC` environment variable to the name or path of
a specific MPI C compiler, if you have many to choose from:

    $ env MPICC="/usr/local/bin/mpicc" ./configure

About support for the GNU Readline library:
------------------------------------------------------------------------

MrBayes will be compiled with support for the [GNU Readline
Library](https://cnswww.cns.cwru.edu/php/chet/readline/rltop.html)
(which allows for command line history within an interactive session,
and command completion through pressing the `tab` key) if the `readline`
library and its associated header files are available.  However, it
seems as if the Readline library doesn't work well together with
common MPI implentations ([OpenMPI](https://www.open-mpi.org/) and
[MPICH](https://www.mpich.org/)), so if MrBayes is configured with
support for MPI parallelization, Readline support will be automatically
disabled.

A possible workaround is to use a "Readline wrapper", such as
[rlwrap](https://github.com/hanslub42/rlwrap), which gives you
persistent command line history between sessions, but no command
tab-completion.  To use rlwrap with an MPI-enabled MrBayes, simply start
MrBayes like this:

    $ rlwrap mb


// vim: filetype=markdown
