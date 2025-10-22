# Atmosphere types

An `EarthAtmosphere` object, as the name suggests, characterizes an atmospheric state on Earth with a number of profiles and other quantities to the extent needed by most models to generate radiances.

Convenience functions exist which make the creation of `EarthAtmosphere` objects easier, such as [`create_empty_EarthAtmosphere`](@ref) or [`create_example_atmosphere`](@ref).

The concept of an atmospheric state supported by RetrievalToolbox is one that includes two categories of profiles. One type of profile is associated with the **retrieval grid**, a one-dimensional vertical pressure profile (`pressure_levels`) that dictates e.g. optical property calculations and also determines the pressures at which the gas profile values are defined. The second  type are the **meteorological profiles**, reserved for specific humidity, temperature, altitude and gravity. These are laid out on another pressure grid (`met_pressure_levels`) which can be a completely different one compared to the retrieval grid.

Each profile has its own physical unit with denoted by `_unit`, so e.g. the altitude levels are stored in `.altitude_levels`, and the corresponding unit field is `.altitude_unit`. The type definition restricts the possible units to those that make sense physically, so users can use `u"km"` or `u"m"` as their physical units for the altitude profile, but not e.g. `u"Pa"`, as that is not a valid length-type unit.

When users write a retrieval algorithm that performs retrievals on many scenes, it is generally advised to not create new atmosphere objects as that tends to fill up memory quite fast. Instead, use the [`ingest!`](@ref) function to copy values into the atmosphere object, over-writing existing values. As with most objects in RetrievalToolbox, users must be cautious as it is always possible to change values inside the object.

```@docs
RE.EarthAtmosphere
```