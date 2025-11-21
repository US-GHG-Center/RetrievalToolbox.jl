---
title: "RetrievalToolbox.jl: A Julia package for trace gas retrievals and related applications"
tags:
  - spectroscopy
  - retrieval
  - trace gas
  - greenhouse gas
  - Julia

authors:
  - name: Peter Somkuti
    orcid: 0000-0002-5858-8471
    affiliation: "1,2"

affiliations:
  - name: Earth System Science Interdisciplinary Center, University of Maryland, MD, USA
    index: 1
  - name: NASA Goddard Space Flight Center, MD, USA
    index: 2

date: 18 November 2025
bibliography: paper.bib
---

# Summary

[`RetrievalToolbox.jl`](https://github.com/US-GHG-Center/RetrievalToolbox.jl) is a software library written in Julia which allows users to write their own trace gas retrieval algorithms and related applications. It provides many types with an established type hierarchy, with each type representing a tangible concept often used in this category of remote sensing problems. Further, `RetrievalToolbox` ships with many functions that act on or utilize those types and provide means to perform various calculations that are common in the field.

A distinguishing feature of `RetrievalToolbox` is that it does not specifically implement a so-called *forward model* (the computation of synthetic observations), giving users the freedom to create their own and be in full control of the main body of calculations.

# Statement of need, related software and alternatives

Passive remote sensing of trace gas concentrations (such as carbon dioxide or methane) in Earth's atmosphere requires some **retrieval algorithm** to extract the desired information from spectroscopic measurements of (reflected) sunlight. Most software implementations of these algorithms are highly bespoke solutions that only work for a particular instrument or set of input data. Despite many of the core assumptions being very similar across various algorithms, adapting an existing implementation to changing user needs proves to be difficult. This is usually due to following reasons: (1) algorithms are designed for performance and a narrow scope of operations, making them hard to read and modify, (2) re-usability and extensive documentation are not funded priorities, and, (3) a high entry barrier for researchers that are less familiar with more complex software architecture.

While most retrieval software that implement specific algorithms are still closed-source, some have been made public. The operational retrieval algorithm for NASA's OCO-2 [@Crisp2004;@Crisp2014] and OCO-3 [@Eldering2019] missions, named ACOS [@Odell2018], is released as [`RTRetrievalFramework`](https://github.com/nasa/RtRetrievalFramework/) -- a highly customizable and modular framework which combines C++ and Fortran routines that are interwoven through Lua scripts. With `RTRetrievalFramework`, users can adjust the behavior of the program flow via the Lua-based umbrella, without having to re-compile the large codebase. A more recent development is the [`ReFRACtor`](https://github.com/ReFRACtor) software, which is a fork of `RtRetrievalFramwork` and replaces the Lua portions with Python scripts.

Hosted at the Netherlands Institute for Space Research (SRON), [`RemoTeC`](https://bitbucket.org/sron_earth/remotec_general/src/main/) [@Butz2011] is the scientific algorithm that forms the basis of the operational methane data products derived from ESA's Sentinel-5P mission [@Veefkind2012]. It is a more traditional software written in Fortran 90 that is controlled at runtime by text-based configuration files. In `RemoTeC`, the forward model is a static part of the architecture of the software, in which the configuration file contents determine which features are active during execution.

The [`QDOAS`](https://github.com/UVVIS-BIRA-IASB/qdoas) application suite focuses on a specific type of gas retrievals (differential optical absorption spectroscopy) and provides both a graphical user interface as well as command line tools for batch processing of measurements. This software is highly popular for use with the MAX-DOAS type of ground-based instruments and provides a streamlined out-of-the-box experience with a helpful user guide. While offering a large number of features, changing the core behavior of the `QDOAS` model requires editing the C/C++ source itself and propagating those changes into the GUI is only possible with knowledge of the [`Qt`](https://www.qt.io/product/framework) framework.

Above examples highlight cases of existing software that cover different application scenarios. They share that newcomers may find it prohibitively difficult to adapt them to their specific use-cases if adjustments are needed. With only highly specialized and often closed-source tools available, it has been challenging to recruit new talent or have already established researchers implement new ideas.

# Design intentions

`RetrievalToolbox` aims to fix the above mentioned obstacles via following approach.

First, it provides object types and functions that can be freely arranged into a retrieval algorithm. It is part of the core design of this software that the explicit algorithm implementation is left up to the user, which is a distinguishing feature compared to existing solutions. This allows users a large degree of flexibility in deciding how the various components interact or which components to include in the first place. For example, when working with an instrument with no sensitivity to the state of polarization, there is no need to represent polarization components that will never be utilized. This flexibility has to be balanced carefully against against overall legibility of the resulting user code. Hiding a series of lengthy expressions in functions may result in a cleaner look for the top-level user code, but will likely also make it more difficult to make non-trivial modifications. `RetrievalToolbox` uses a balanced approach. Some tasks, such as calculating altitude and gravity profiles can be done in a very compact manner

``` julia
# Calculate the altitude and gravity profiles of `scene`,
# an `EarthScene` object.
calculate_altitude_and_gravity!(scene)
```

whereas others are more verbose, like copying over the freshly computed model radiances into the expected place, while accounting for the correct units:

``` julia
# Copy result from inst_buf -> rt_buf
@views rt_buf.radiance.I[rt_buf.indices[swin]] =
    inst_buf.low_res_output[dispersion[spec].index] *
    buf.rt[swin].radiance_unit / rt_buf.radiance_unit
```

More verbose expressions can assist in making certain aspects of a user program more clear, but should be used sparingly.

Julia, as the programming language of choice [@Bezanson2017] in which `RetrievalToolbox` was implemented, plays an important role. With Julia it is relatively straightforward to write procedures in a high-level scripting manner, while still retaining the performance of compiled code. Further, Julia programs or smaller program snippets are easily run in notebooks such as Jupyter [@Kluyver2016] or Pluto [@Fons_van_der_plas_2025], which make authoring interactive training material straightforward.

Second, the source code itself is extensively documented. Functions and types contain documentation via *docstrings*, which provide both helpful context regarding the motivation and intended usage as well as the explicit function interface. Docstrings can be accessed interactively through the Julia REPL. For example, typing

``` julia
julia> using RetrievalToolbox
julia> ?create_pressure_weights
```

will bring up helpful information regarding the `create_pressure_weights` function. In addition to dosctrings, the code repository also hosts online documentation in the form of a website. Built with [`Documenter.jl`](https://documenter.juliadocs.org/), the documentation lists the types and functions along with their docstrings grouped into appropriate categories. It includes details on overall usage, design principles and code snippets, as well as guides on how to extend the software library with new types and functions if the provided ones do not fully meet user needs. Additionally, several example implementations already exist and are linked on the [`Github page`](https://github.com/US-GHG-Center/RetrievalToolbox.jl).

Finally, users can read through [`learning materials`](https://petersomkuti.github.io/RetrievalToolbox-Tutorials/) at their own pace. These tutorials provide an entry point for researchers with the relevant background for learning how to use `RetrievalToolbox`.

# Scope and functionality

`RetrievalToolbox` implements concepts that are common in many areas of passive remote sensing. At the core is a vertically layered model atmosphere in which atmospheric constituents (gases, aerosols, ..) are considered horizontally homogeneous. Relevant total optical properties are calculated based on the particular objects present in the model atmosphere, which are then ingested by the radiative transfer functions to produce model radiances. For absorption-only set-ups, internal routines are available to perform radiative transfer calculations based on the Beer-Lambert-Bouguer law. If atmospheric scattering (aerosols, Rayleigh molecular scattering) must be accounted for, `RetrievalToolbox` supports interfacing with the [`XRTM`](https://github.com/gmcgarragh/xrtm) radiative transfer library.

To facilitate the creation of retrieval algorithms, `RetrievalToolbox` contains a number of pre-defined state vector element types (e.g. surface reflectance polynomial coefficients, gas volume mixing ratios, etc.) which trigger internal calculation of partial derivatives. Inverse solver constructs then use those derivatives to compute Jacobian matrices that are required to iteratively adjust the state vector elements in order to match some measurement.

Below follows a short (incomplete) list of some currently supported features

* Single-scene paradigm with 1-dimensional radiative transfer through a layered model atmosphere
* Suitable for the wavelength range between UV (~200 nm) and upper mid-IR (~2700 nm)
* Dedicated types to represent scalar-only or vector-valued radiances (I,Q,U components)
* Arbitrary spectral window configuration (single-band, multi-band)
* Down- and up-looking instrument viewing geometries
* Built with optimal-estimation type [@Rodgers2000inverse] inversions in mind
* Uses [`Unitful.jl`](https://github.com/JuliaPhysics/Unitful.jl) to represent physical units

# Acknowledgements

`RetrievalToolbox` received additional funding through NASA's Cooperative Earth Systems Science Research Agreement (CESSRA), grant number 80NSSC23M0011.

# References

