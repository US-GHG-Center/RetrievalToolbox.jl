# [Buffer types] (@id buffer_types)

Below are the currently supported buffer types which represent collections of pre-allocated objects that can be used over and over in situations where users want to e.g. perform many retrievals. The following buffers can be created once and then used repeatedly for different scenes. The slower the retrieval process itself is, the lower is the overall performance gain from re-using buffers, however. For very fast retrievals with completion times less then a few seconds, it is highly recommended to re-use the buffers.

## How to construct and use an `EarthAtmosphereBuffer`

The `EarthAtmosphereBuffer` type requires itself several other objects based on RetrievalToolbox types. Most importantly, a number of functions that RetrievalToolbox provides implicitly assume that all those objects are correctly instantiated and all relevant relationships hold across the type hierarchy. While it is possible to create an `EarthAtmosphereBuffer` object manually, it is highly recommended to use the convenience function that is provided with RetrievalToolbox. In cases where this convenience function does not produce the right configuration due to specific user needs, the best course of action is to copy the function code (`src/buffers.jl`) into a new module and make appropriate changes to capture those cases. See also: [How to develop or extend RetrievalToolbox](@ref develop_RetrievalToolbox).

Before an `EarthAtmosphereBuffer` can be successfully created, users must first have various other objects prepared beforehand.

* `sv`: some `AbstractRetrievalStateVector` objects, which can be empty, i.e. write `sv = RetrievalStateVector([])`. Note the following, however.
* `spectral_windows`: a list of `AbstractSpectralWindow` objects
* `surface_types`: a list of tuples, indicating the type of surface to be used, along with parameters (see below: [surface types](@ref surface_types_buffer))
* `atmospheric_elements`: a list of `AbstractAtmosphereElement` objects
* `solar_models`: a dictionary which maps `AbstractSpectralWindow` objects to `AbstractSolarModel` ones so that RetrievalToolbox knows which solar model is to be used for which spectral window. Several (or all) spectral windows may point to the same solar model.
*  `RT_models`: a list of symbols to indicate which type of radiative transfer model is to be used for each spectral window. The order **must be the same** as given in the `spectral_windows` list.
*  `RadType`: the type of radiance to be used, either `ScalarRadiance` or `VectorRadiance`
*  `rt_buf`: an `AbstractRTBuffer`
*  `inst_buf`: an `InstrumentBuffer`
*  `N_level`: the number of retrieval levels of the [`EarthAtmosphere`](@ref)
*  `N_met_level`: the number of meteorology profile levels of the [`EarthAtmosphere`](@ref)
*  `T`: the number type to be used for arrays (`Float64` is recommended)


```@docs
RE.EarthAtmosphereBuffer(
    sv::AbstractStateVector,
    spectral_windows,
    surface_types::Vector{Tuple},
    atmospheric_elements,
    solar_models::Dict{<:AbstractSpectralWindow,<:AbstractSolarModel},
    RT_models::Vector{Symbol},
    RadType::Type{<:Radiance},
    rt_buf::AbstractRTBuffer,
    inst_buf,
    N_level::Integer,
    N_met_level::Integer,
    T::Type{<:AbstractFloat}
    )
```

!!! warning
    In **most** cases, users will want to make sure that `sv`, the state vector supplied to this convenience function, is the same state vector that is used to map the keys of the `jacobian` dictionary of the RT buffer `rt_buf`. There is currently **no check in place to verify this**.


### [Surface types](@id surface_types_buffer)

Both `RT_models` and `surface_types` arguments must be ordered according to the `spectral_windows` argument (list of spectral windows). For example, let there be three spectral windows we want to use, `swin_A`, `swin_B`, `swin_C` such that `spectral_windows = [swin_A, swin_B, swin_C]`. Users have the choice to perform RT calculations with different models for each spectral window. So if the first two windows should be run with a Beer-Lambert-type RT and the last one use the XRTM library, one would write `RT_models = [:BeerLambert, :BeerLambert, :XRTM]`. Similarly, we proceed to assign surface types. At the moment only Lambertian-type surface are supported with the [`BeerLambertRTMethod`](@ref), but for a [`MonochromaticRTMethod`](@ref) that employs the XRTM solver, we can make use of the RPV BRDF kernel:

```julia
surface_types = [
    (:Lambertian, 2),
    (:Lambertian, 2),
    (:RPV, 2, 0.05, -0.1, 0.75),
    ]
```

For the `Lambertian` surface, there is only one needed parameter, which is the order of the polynomial to capture the spectral dependence. For the `RPV` BRDF kernel, three additional parameters need to be part of the tuple: the hotspot, asymmetry and anisotropy parameters (see also [`RPVPolynomialKernel`](@ref)). The important part here is that the surface are created and assigned, in order such that `swin_A` and `swin_B` both get their own Lambertian surfaces, and the RPV BRDF kernel surface will be attached to `swin_C`.



## Types

```@docs
RE.EarthAtmosphereBuffer
```

```@docs
RE.InstrumentBuffer
```

```@docs
RE.ScalarRTBuffer
```

```@docs
RE.VectorRTBuffer
```