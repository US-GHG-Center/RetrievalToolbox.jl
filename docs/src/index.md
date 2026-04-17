RetrievalToolbox is a software library written in pure Julia to facilitate building trace gas retrieval algorithms and related applications.

As a first stop, reading through the [Fundamentals](@ref) section will provide both a broad overview of the usage principles of RetrievalToolbox.

The [Design](@ref) section lays out the design philosophy of the software and can help new users to understand specific workings of RetrievalToolbox or how different functions interact with each other. Users will find explanations on how we utilize abstract types, how to deal with wavenumbers and wavelengths, how physical units are handled (and where not), and how users should think about writing a retrieval algorithm with RetrievalToolbox.

To users that are (relatively) new to Julia, we have also integrated a few documents that explains various concepts of Julia that are often used in RetrievalToolbox, such as the usage of dictionaries (see [here](@ref working_with_dicts)), which may differ from other programming and scripting languages, or how physical units are considered (see [here](@ref working_with_units)).


In the current version of the software library, the following are not (yet) supported:

* Spectrally varying aerosol scattering properties:
  * At the moment, aerosol scattering properties, as represented by the phase matrix (or phase function) expansion coefficients, are evaluated only for the center of spectral bands.
* Limb soundings

!!! warning
    API documentation is not complete at this time. We recommend users to study the example implementations found on the main repository page, reach out directly or via the [GitHub Issues Tracker](https://github.com/RetrievalToolbox/ACOS-Goddard/issues).