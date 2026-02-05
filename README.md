# RetrievalToolbox
[![docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Static Badge](https://img.shields.io/badge/Tutorials-green)](https://retrievaltoolbox.github.io/RetrievalToolbox-Tutorials/)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

RetrievalToolbox is a library for building trace gas retrieval algorithms and related applications written in pure [Julia](https://julialang.org). The library is currently in an early release stage, and feature-breaking updates might happen - although we attempt to keep those to a minimum. For the time being, we recommend to fork this repository into your own GitHub organization and integrate updates from here as to not break your own application.

RetrievalToolbox was developed at the Earth System Science Interdisciplinary Center (ESSIC) at the University of Maryland College Park, and at NASA Goddard Space Flight Center.


## Documentation and Learning

The main documentation, built with [Documenter.jl](https://documenter.juliadocs.org/stable/) is part of the repository under `docs/`, and the corresponding HTML render can be found [here][docs-dev-url]. Learning materials can be found [here](https://retrievaltoolbox.github.io/RetrievalToolbox-Tutorials/) - new users are **strongly encouraged** to read through these tutorials.

Looking at working examples is also highly instructive:

 * [EMIT retrieval demo for AGU2025 (New Orleans, USA, 2025)](https://github.com/RetrievalToolbox/EMIT-retrieval/)
 * [Demo for IWGGMS21 (Takamatsu, Japan, 2025)](https://github.com/US-GHG-Center/IWGGMS21-Demo)
 * [Implementation of NASA's ACOS algorithm](https://github.com/RetrievalToolbox/ACOS-Goddard/)

_Users are very welcome to submit their working set-ups to be listed here!_

## Installation

The library can be installed directly from Julia by typing (via http)

    using Pkg
    Pkg.add(url="https://github.com/US-GHG-Center/RetrievalToolbox.jl")

or (via ssh)

    using Pkg
    Pkg.add(url="git@github.com:US-GHG-Center/RetrievalToolbox.jl")

This will install the RetrievalToolbox including all needed dependencies.

If you forked this repository, you must amend the above commands to pull the package from your new repository location.


### Building XRTM and making RetrievalToolbox aware of its location

RetrievalToolbox makes use of the XRTM library to perform the various radiative transfer calculations which require scattering from, e.g. molecular Rayleigh scattering or aerosols. While it is possible to use RetrievalToolbox with a built-in Beer-Lambert-Bouguer method, many retrieval applications will need to account for scattering and thus require the XRTM library. XRTM is published at (https://github.com/gmcgarragh/xrtm/), and we also maintain a fork at (https://github.com/PeterSomkuti/xrtm). Depending on your needs and the way how you create your own algorithm, you might prefer one vs. the other option.

To build XRTM, you need somewhat recent versions of GCC, gfortran and make. First, clone the repository into some location on your computing environment:

`git clone https://github.com/gmcgarragh/xrtm`

or

`git clone https://github.com/PeterSomkuti/xrtm`

Note that as of May 2025, the first repository does not compile on MacOS without modifications to the makefile. The second repository, however, should compile successfully. In either case, you must take a copy of the example makefile and save it as `make.inc`. So switch to the XRTM directory (`cd xrtm`) and

`cp make.inc.example make.inc`

Now you have to edit `make.inc` to work seamlessly with your environment. Usually, you will need at least the following edits:

1. `CC= ` set your (GNU) C compiler executable. Something like `CC=gcc` might work on most computers.
2. `CXX= ` set your (GNU) C++ compiler executable. `CC=g++` should work.
3. `F77= ` and `F90= ` set your (GNU) Fortran77 and Fortran90 compilers. `F77=gfortran` and `F90=gfortran` should work.
4. Set your LAPACK library location. On MacOS, this likely needs to be set to `-framework Accelerate`, and on Linux machines, you will likely have many options, including various OpenBLAS options, or Intel's MKL.

**On recent versions of MacOS, `gcc` invokes the Clang compiler rather than the GNU compiler suite.** Apple Clang does not work with the `-fopenmp` flag, and thus building XRTM will not work. We recommend downloading GCC (and gfortran) via [`homebrew`](https://brew.sh), by typing `brew install gcc`. To then call this freshly installed GCC compiler, one needs to add the version number, i.e. `gcc-15`. Use `gcc-15`, `g++-15` etc. to edit the Makefile.

Once those are set, simply type `make` to build the XRTM library. Note that you will have to make more substantial changes if you want to compile with e.g. Intel's compiler suite. If you do not have a LaTeX distribution with `pdflatex` on the computer, the build process will exit with an error message, but that only concerns the compilation of the documentation which is part of the makefile; by that point the library will have been built successfully and is ready to use.

RetrievalToolbox looks in a particular place to find the XRTM Julia interface file (which is part of XRTM), and it is determined by the environment variable `XRTM_PATH`, which should point to the main source tree that you just cloned. So in an interactive environment, you should, for example, start a Julia session with

`XRTM_PATH=/path/to/xrtm julia`

Or if you use RetrievalToolbox in a jupyter environment, start it up appropriately with e.g.

`XRTM_PATH=/path/to/xrtm jupyter-lab`

Alternatively, if that variable is exported beforehand, the variable will be also correctly set in Julia. Environment variables in Julia are accessible in the `ENV` dictionary, which is also writable from within a Julia session. So if `XRTM_PATH` is not set before the Julia session started, you can simply change it before loading the RetrievalToolbox module:

``` julia
ENV["XRTM_PATH"] = "/path/to/xrtm"
using RetrievalToolbox
```

**Important!** The second repository (https://github.com/PeterSomkuti/xrtm) contains a small, but very impactful code modification. In `xrtm/src/xrtm.h` the pre-processor directive was changed to say
``` c
#define DO_NOT_ADD_SFI_SS 1
```
(instead of `#define DO_NOT_ADD_SFI_SS 0`). This has the effect that for certain RT solvers, such as the `eig_bvp` or `two_stream` ones, the contributions from single-scattering **are not automatically computed and added** to the total radiance fields and their derivatives. This is the wanted behavior for some applications, such as retrievals from NASA's OCO instruments, where you may want to compute the single-scatter contributions with a vector RT call, but the diffuse (MS) contributions with a scalar RT call. **Be aware that this is a compile-time choice** at the moment, so switching between `#define DO_NOT_ADD_SFI_SS 0` and `#define DO_NOT_ADD_SFI_SS 1` requires re-compiling. Alternatively, you can keep two copies of the code with the two different variants for this variable, and point RetrievalToolbox to a different path when you run it via changing `XRTM_PATH`.

## Citing RetrievalToolbox

We will be working in the future to submit an entry to the Journal of Open Source Software (JOSS), which will allow citing RetrievalToolbox with a proper DOI.

## Contributing

Contributions are always welcome, whether bug reports, bug fixes (through pull requests), suggestions, feature requests or otherwise. Please use the Github [issue tracker](https://github.com/RetrievalToolbox/ACOS-Goddard/issues).

## Alternatives

RetrievalToolbox was heavily influenced by the hard work of numerous scientists in different labs and institutions in various countries. Some of those algorithms have been under development for over a decade and a proven track record of reliability. Prominent and publicly available alternatives are

- NASA's RtRetrievalFramework (https://github.com/NASA/RtRetrievalFramework) and ReFRACtor (https://github.com/ReFRACtor)
- SRON's RemoteC (https://bitbucket.org/sron_earth/remotec_general/src/main/)


[docs-dev-url]: https://US-GHG-Center.github.io/RetrievalToolbox.jl/dev/
