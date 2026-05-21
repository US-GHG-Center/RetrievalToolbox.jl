# [Aerosol types](@id aerosol_types)

Aerosols are part of the `AbstractAtmosphereElement` type tree: `AbstractAtmosphereElement` → `AbstractAerosolType`, and have to be added to the model atmosphere into the `atm_elements` vector of an [`EarthAtmosphere`](@ref). Some aerosol type belonging to `AbstractAerosolType` is a more specific type implementation of an aerosol type. Aerosols may need optical properties, which a grouped into the `AbstractAerosolProperty` abstract types.

Generally, an `AbstractAerosolType` represents the abundance and vertical distribution of an aerosol (or aerosol mixture), and the attached `AbstractAerosolProperty` will contain the information about the optical properties such as the aerosol's scattering characteristics. Aerosol types are mutable structures, so users can change their values on-the-fly or through state vector elements.

```@docs
RE.GaussAerosol
```

```@docs
RE.MieMomAerosolProperty
```