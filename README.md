
RetrievalToolbox is a library for building trace gas retrieval algorithms and related applications written in pure [Julia](https://julialang.org).

# Installation

The library can be installed directly from Julia by typing (via http)

    using Pkg
    Pkg.add(url="https://github.com/US-GHG-Center/RetrievalToolbox.jl")

or (via ssh)

    using Pkg
    Pkg.add(url="git@github.com:US-GHG-Center/RetrievalToolbox.jl")

This will install the RetrievalToolbox including all needed dependencies.

## Building XRTM and making RetrievalToolbox aware of its location

RetrievalToolbox makes use of the XRTM library to perform the various radiative transfer calculations which require scattering from, e.g. molecular Rayleigh scattering or aerosols. While it is possible to use RetrievalToolbox with a built-in Beer-Lambert-Bouguer method, many retrieval applications will need to account for scattering and thus require the XRTM library. XRTM is published at (https://github.com/gmcgarragh/xrtm/), and we also maintain a fork at (https://github.com/PeterSomkuti/xrtm). Depending on your needs and the way how you create your own algorithm, you might prefer one vs. the other option.

To build XRTM, you need somewhat recent versions of GCC, gfortran and make. First, clone the repository into some location on your computing environment:

`git clone https://github.com/gmcgarragh/xrtm`

or

`git clone https://github.com/PeterSomkuti/xrtm`

Note that as of May 2025, the first repository does not compile on MacOS without modifications to the makefile, the second repository, however, should compile successfully. In either case, you must take a copy of the example makefile and save it as `make.inc`. So switch to the XRTM directory (`cd xrtm`) and

`cp make.inc.example make.inc`

Now you have to edit `make.inc` to work seamlessly with your environment. Usually, you will need at least the following edits:

1. `CC= ` set your (GNU) C compiler executable. Something like `CC=gcc` might work on most computers.
2. `CXX= ` set your (GNU) C++ compiler executable. `CC=g++` should work.
3. `F77= ` and `F90= ` set your (GNU) Fortran77 and Fortran90 compilers. `F77=gfortran` and `F90=gfortran` should work.
4. Set your LAPACK library location. On MacOS, this likely needs to be set to `-framework Accelerate`, and on Linux machines, you will likely have many options, including various OpenBLAS options, or Intel's MKL.

Once those are set, simply type `make` to build the XRTM library. Note that you will have to make more substantial changes if you want to compile with e.g. Intel's compiler suite.

RetrievalToolbox looks in a particular place to find the XRTM Julia interface file, and it is determined by the environment variable `XRTM_PATH`, which should point to the main source tree that you just cloned. So in an interactive environment, you should, for example, start a Julia session with

`XRTM_PATH=/path/to/xrtm julia`

Alternatively, if that variable is exported beforehand, the variable will be also correctly set in Julia. Environment variables in Julia are accessible in the `ENV` dictionary, which is also writable from within a Julia session. So if `XRTM_PATH` is not set before the Julia session started, you can simply change it before loading the RetrievalToolbox module:

``` julia
ENV["XRTM_PATH] = "/path/to/xrtm"
using RetrievalToolbox
```

## Learning

We are working on finalizing a set of tutorials for new users that introduce the basic concepts of RetrievalToolbox. Once ready, they will be linked here.

## Citing RetrievalToolbox

We will be working in the future to submit an entry to the Journal of Open Source Software (JOSS), which will allow a citing RetrievalToolbox with a proper DOI.

## Alternatives

RetrievalToolbox was heavily influenced by the hard work of numerous scientists in different labs and institutions in various countries. Some of those algorithms have been under development for over a decade and a proven track record of reliability. Prominent and publicly available alternatives are

- NASA's RtRetrievalFramework (https://github.com/NASA/RtRetrievalFramework)
- SRON's RemoteC (https://bitbucket.org/sron_earth/remotec_general/src/main/)