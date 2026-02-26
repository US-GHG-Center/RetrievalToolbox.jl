# [Surface types](@id surface_types)

Surface types are required by the scene objects (`EarthScene`) to carry information about what type of surface model the radiative transfer calculations should use. As laid out in the `EarthScene` type, one surface is mapped to exactly one spectral window - there is currently no support for "sharing" a surface between different spectral windows. If users want to fix the retrieved surface properties between two separate retrieval windows, they must manually impose correlations through the prior covariance matrix.

The [`NoSurface`](@ref) type can be used for up-looking instrument retrievals (e.g. sun-tracking spectrometers).

The [`BeerLambertRTMethod`](@ref) currently supports only the [`NoSurface`](@ref) and [`LambertianPolynomialSurface`](@ref) surface types. The [`MonochromaticRTMethod`](@ref) using the XRTM solver requires a different set of surface properties as they are set-up to be potentially polarized BRDFs with strong angular dependence, and they can be recognized by the "Kernel" suffix in the surface type name. For example, [`LambertianPolynomialKernel`](@ref) represents a (spectrally varying) Lambertian surface in the XRTM context, and [`RPVPolynomialKernel`] is a (spectrally varying) surface characterized by the Rahman-Pinty-Verstraete[^RPV1993] kernel.



```@autodocs
Modules = [RetrievalToolbox]
Pages = ["surface_types.jl"]
Order = [:type]
```

[^RPV1993]: [Rahman, H., B. Pinty, and M. M. Verstraete (1993), Coupled surface-atmosphere reflectance (CSAR) model: 2. Semiempirical surface model usable with NOAA advanced very high resolution radiometer data, J. Geophys. Res., 98(D11), 20791–20801](https://doi.org/10.1029/93JD02072)