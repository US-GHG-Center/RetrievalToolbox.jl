# Atmosphere types

An `EarthAtmosphere` object, as the name suggests, characterizes an atmospheric state on Earth with a number of profiles and other quantities to the extent needed by most models to generate radiances.

Convenience functions exist which make the creation of `EarthAtmosphere` objects easier, such as [`create_empty_EarthAtmosphere`](@ref) or [`create_example_atmosphere`](@ref).

The concept of an atmospheric state supported by RetrievalToolbox is one that includes two categories of profiles. One type of profile is associated with the **retrieval grid**, a one-dimensional vertical pressure profile (`.pressure_levels`) that dictates e.g. optical property calculations and also determines the pressures at which the gas profile values are defined. The second  type are the **meteorological profiles**, reserved for specific humidity, temperature, altitude and gravity. These are laid out on another pressure grid (`.met_pressure`) which can be a completely different one compared to the retrieval grid.

Note that retrieval grid explicitly differentiates `pressure_levels` and `pressure_layers`, whereas the meteorological grid does not. This is done to reduce potential confusion about the double meaning of the term "level". For the purposes of radiative transfer and retrievals, we use "level" explicitly to refer to a layer boundary, where "layer" is indeed a horizontal slab with a finite thickness. Gas volume mixing ratios are defined on levels, rather than layers, as it makes radiative transfer and Jacobian calculations easier.

The meteorological profiles (temperature, specific humidity, altitude, gravity) are usually taken from NWP models such as NASA's GEOS or ECMWF's IFS systems, where the relevant quantities are defined on **model levels**. These model levels, however are often horizontal slabs with finite thickness (as opposed to our notion of level being a layer boundary), with a reported mid-level pressure.

!!! tip
    We highly recommend to ingest meteorological profiles from NWP model outputs **as-is**, by setting the `.met_pressure` grid equal to the model mid-level pressures, and then simply copying the relevant profiles (temperature → `.temperature`, specific humidity → `.specific_humidity`, ..). RetrievalToolbox internally uses interpolation to obtain the profile values at arbitrary pressures when needed for various calculations (e.g. optical depth integration). This provides smooth profiles, as long as the retrieval grid is roughly similar or coarser than the meteorological grid.


Each profile has its own physical unit with denoted by `_unit`, so e.g. the altitude profile are stored in `.altitude`, and the corresponding unit field is `.altitude_unit`. The type definition restricts the possible units to those that make sense physically, so users can use `u"km"` or `u"m"` as their physical units for the altitude profile, but not e.g. `u"Pa"`, as that is not a valid length-type unit.

When users write a retrieval algorithm that performs retrievals on many scenes, it is generally advised to not create new atmosphere objects as that tends to fill up memory quite fast. Instead, use the [`ingest!`](@ref) function to copy values into the atmosphere object, over-writing existing values. As with most objects in RetrievalToolbox, users must be cautious as it is always possible to change values inside the object.

```@docs
RE.EarthAtmosphere
```