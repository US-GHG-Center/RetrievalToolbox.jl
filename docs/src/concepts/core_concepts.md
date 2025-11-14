# [Core Algorithm Concepts](@id core_concepts)

!!! note
    Within the RetrievalToolbox algorithm tools, an attempt is made to keep notation consistent as well as stick to more modern terminology. This approach occasionally overrides more antiquated terms, so various sections in the documentation will emphasize when a possible clash is expected.

## Instrument Spectral Response Function (formerly ILS)


[OCO ATBD](https://docs.julialang.org/en/v1/stdlib/Markdown/)

## Pixels, Spectral Samples and Dispersion

The distinction between pixels, samples and what dispersion describes occasionally causes confusion due to a lack of consistent terminology in publications or other various documents.

A pixel is considered to be a discrete unit on an instrument detector of any type (where appropriate). Some instruments aggregate pixels during the read-out process such that the data received does not truly reflect the physical detector elements. Other instruments do read out the detector on a native pixel level, however some form of processing is performed afterwards and the resulting data, as ingested by retrieval algorithms, can no longer considered to be per-pixel.

Retrieval algorithms, in general, act on calibrated, geo-located radiance data, often denoted as Level-1b (L1B, L1b). For hyperspectral data, there usually is at least one spectral dimension of that data such that when the data is extracted along that dimension, one obtains what is considered a **spectrum**. One spectrum thus has a discrete number of elements along its spectral axis and the spectral axis consists of spectral samples. Those spectral samples do not have to be identical with the underlying detector pixels. For most instruments, they are not. Even for instruments in which the detector pixels map 1:1 into spectral samples, it is a good choice to stay consistent in the terminology and refer to an element of a spectrum as spectral sample.

For a single spectrum ``I`` that is extracted from calibrated, geo-located radiance data, the spectral sample information can be written explicitly as ``I_{[s]}``, with $s$ being some index of its spectral dimension. ``I_{[s]}`` is naturally a discrete quantity with $s$ being a discrete index itself. Compare this to a theoretical description of radiance, which in general is a continuous function of wavelength or wavenumber ``\tilde{I}(\lambda)`` or ``\tilde{I}(\nu)``. In order to allow for a comparison between measured quantity ``I_{[s]}`` and a model radiance, one must know which wavelength or wavenumber corresponds to a spectral sample at index ``s``.

The relationship between spectral sample (index) and wavelength or wavenumber is usually called **dispersion**. An alternative term for dispersion is **wavelength (wavenumber) solution**. It is some general function $d$ which maps a spectral sample index to either wavelength or wavenumber, whichever is appropriate for the specific instrument: ``d(s) = \lambda_{[s]}`` or ``d(s) = \nu_{[s]}``. When ``d`` is known, it is straightforward to evaluate some continuous function ``\tilde{I}`` at the correct wavelength or wavenumber in order to compare it to a measured value: ``I_{[s]} \sim \tilde{I}{\left( d(s) \right)}``.

For many instruments, the function ``d`` is generally smooth and tends to be expressed as a polynomial which maps spectral sample to wavelength or wavenumber in the following way:

```math
    \lambda_{[s]} = \sum_{i=0}^{N} c_i \cdot s^i, \\
```

or

```math
    \nu_{[s]} = \sum_{i=0}^{N} c_i \cdot s^i.
```

The polynomial coefficients ``c`` are usually either available in the published measurement data, in accompanying documents, or in rare cases, have to be derived.