# State Vector Types

Certain components of RetrievalToolbox require an `AbstractStateVector` that controls the computation of quantities relevant to a retrieval algorithm. For example, during the calculation of optical properties, the routine scans through the supplied `AbstractStateVector` to decide whether certain partial derivatives are needed.

State vector types can be understood as custom structures that only contain the most necessary information about them as needed in an optimal estimation-type retrieval algorithm, and variables of these types are self-descriptive. State vector elements dictate the behavior of functions that are employed in forward models. For example, the `apply_radiance_correction!` function is implemented for appropriate state vector elements. Somewhere in the forward model, it will be appropriate to alter the forward model calculation, depending on the state vector.

The code section below shows a simple example in which some radiance correction is applied to a state vector object `state_vector`. In this example, a loop runs through all state vector elements and the appropriate `apply_radiance_correction!` function is dispatched on the specific type of the state vector element object.

```julia
# My forward model

# [...]
# Calculate radiances here
# [...]

# Apply radiance correction, depending on state vector
for sve in state_vector.state_vector_elements
    RE.apply_radiance_correction!(buffer, sve)
end

# [...]

```

This way of writing a forward model makes use of Julia's multiple dispatch paradigm via which the most specific function will be called, if it is implemented. The strong benefit to users is that the forward model itself remains very legible, as no long conditional statements are needed to determine when to perform the required computation. The default implementation of the function in this example does nothing, so any state vector elements that should not induce any change in the forward model at this point, will act as expected.

More fundamentally, any forward model must incorporate the state vector its calculations through updating the state upon which the forward model acts. **Users must implement necessary functionality themselves.**  Since RetrievalToolbox does not provide a forward model per se, but rather provides the many components that users can utilize to build one, there are many functions which require a state vector type, as listed further in [State Vector Functions](@ref).

## Retrieval State Vector

A `RetrievalStateVector` can be instantiated by passing a vector of  `AbstractStateVectorElement` variables.

```@docs
RE.RetrievalStateVector
```

# State Vector Element (SVE) Types

An `AbstractStateVector` consists of a list of state vector element (SVE) types, which are all subtypes of `AbstractStateVectorElement`. Every `AbstractStateVectorElement` type must be implemented as a `mutable struct`, due to the overall design of RE. Since instantiating new state vectors and state vector elements for each new retrieval scene would impact highly negatively on the overall performance, mutable structs are used. Once instantiated, users must make sure to reset each state vector element at the beginning of a new retrieval by setting the appropriate values and emptying out the `iterations` vector field with `empty!(sve.iterations)`, with `sve` being the state vector element variable.

## Gas scaling factors

For certain retrieval applications, it is sufficient to scale the volume mixing ratio profile of some target gas, rather than having to retrieve the full shape. The `GasLevelScalingFactorSVE` allows to define a section of the volume mixing ratio profile (or the full profile from top to bottom), which is adjusted during the retrieval.

If, for example, in a 10-level atmosphere set-up `start_level` and `end_level` are chosen to be (1, 5), then this SVE will only change the upper half of the gas VMR profile. To let this SVE scale the full profile, in this case the two values should be set to (1, 10), as the values are inclusive - inclusive meaning that the scale factor affects VMR profiles at levels 1,2,3,..,8,9,10. Multiple `GasLevelScalingFactorSVE` elements can appear in a state vector, each referencing a different portion of the atmosphere. As of now, there is no error catching mechanism, so users could define overlapping sections that will likely result in poor retrieval performance.

Within the calculations of optical properties, there are checks in place to determine whether partial derivatives of the optical properties with respect to volume mixing ratios are required. If a `GasLevelScalingFactorSVE` is part of the state vector, the partial derivatives of optical depth with respect to volume mixing ratio is automatically calculated.

This type requires a unit without dimensions, such as `Unitful.percent` or `Unitful.NoUnits`.

!!! note
    In order for gas objects to be correctly updated, users must call the `atmosphere_element_statevector_update!` function as part of their forward model. Further, it is important to roll back the updates with the `atmosphere_element_statevector_rollback!` function towards the end of the forward model. For more details, consult following section in the manual: [Updates and roll-backs](@ref).


```@docs
RE.GasLevelScalingFactorSVE
```

## SIF radiance

Solar-induced fluorescence (SIF) can be co-retrieved in some circumstances where the retrieval window covers part of the range of SIF emission, roughly between 670 nm and 850 nm. If the model atmosphere has a `SIFRadiance` object in the list of `AtmosphereElements`, and the `SIFRadianceSVE` points to that object, then Jacobians are calculated and the SIF emission magnitude can be adjusted for during the inversion.

!!! warning
    Note that the radiance units of the `SIFRadianceSVE` state vector element object, and the radiance unit of the `SIFRadiance` object need not be the same, as internal functions make appropriate conversions - however we highly recommend to **not mix** units of Watts and photons / s.

There is also **no specific unit check** for radiance-type units, so users must make sure themselves that the correct radiance-type unit is supplied, otherwise an error will be raised at some point when incompatible units clash during the calculation of the wavelength-dependent SIF radiance.


```@docs
RE.SIFRadianceSVE
```

## Surface pressure

The atmospheric model in RetrievalToolbox, the surface pressure is equal to the lowest retrieval pressure level. Retrieving surface pressure through a `SurfacePressureSVE` allows shifting that lowest level. As with, e.g. the `GasLevelScalingFactorSVE`, the optical property routines will automatically perform the necessary calculations of partial derivatives that are needed to obtain the appropriate Jacobian. This state vector element provides the Jacobian that reflects a positive change of the lowest pressure level, i.e. ``\partial I / \partial p_{N_\mathrm{lev}}``.

Note that there is no default or automatic adjustment of the retrieval pressure level grid that takes into account a possible surface pressure retrieval.

!!! warning
    Users must implement the adjustment of surface pressure into the retrieval grid themselves. A failure to do so will result in poor retrieval performance.

This state vector element requires a unit that is compatible with pressure (e.g. Pascal) which does not need to be the same unit as that of the retrieval pressure level grid. The code internally handles the conversion. Users could, for example, have a retrieval grid in Pascal, but retrieve the surface pressure in units of Torr.

```@docs
RE.SurfacePressureSVE
```

## Surface albedo polynomial

For a Lambertian surface model, in particular a `LambertianPolynomialSurface`, users can
retrieve the polynomial coefficients. This state vector element has to be attached to a spectral window, since surfaces themselves are attached to a spectral window. Like all polynomial-related
quantities, any state vector element must relate to a particular coefficient order. A convenience constructor exists which populates the `unit` field appropriately, given the coefficient order. One might call, for example

```julia
RE.SurfaceAlbedoPolynomial(
    swin, # Spectral window
    2, # Polynomial order
    u"µm", # Associated unit
    0.0, # First guess
    0.0, # Prior value
    0.0, # Prior covariance
)
```
where the appropriate unit `u"µm^-2"` will be generated automatically, and the user must only supply some unit that is compatible with the spectral unit of the spectral window, i.e. a length unit for wavelengths, or a wavenumber unit for wavenumbers. An error will be thrown if units are incompatible.

This state vector element creates a wavelength-dependent quantity and uses the reference wavelength given by the linked spectral window. When calculating the surface reflectivity via `calculate_surface_reflectivity!`, the specialized function `calculate_Lambertian_surface_reflectivity!` computes the per-wavelength albedo using the relative wavelength ``\lambda - \lambda_{\mathrm{ref}}`` or relative wavenumber ``\nu - \nu_{\mathrm{ref}}``.

```@docs
RE.SurfaceAlbedoPolynomialSVE
```

In order to update the surface(s) given some state vector, a helper function exists which performs the in-place update within a `RE.EarthScene` object: `RE.surfaces_statevector_update!`.

## Dispersion polynomial



```@docs
RE.DispersionPolynomialSVE
```

## Zero level offset polynomial

```@docs
RE.ZeroLevelOffsetPolynomialSVE
```