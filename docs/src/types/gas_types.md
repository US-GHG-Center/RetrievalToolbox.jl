# [Gas Types](@id gas_types)

Gases are represented in RetrievalToolbox through their own type `GasAbsorber` which is a sub-type of `AbstractAtmosphereElement`. Each gas must have its an associated spectroscopy table, along with a vector that will contain the volume mixing ratios (VMR) on (pressure) levels, as well as a descriptor.

!!! warning
    RetrievalToolbox routines assume that the VMR level vector is of the same length as the retrieval pressure grid, however no explicit checks are done yet! **You must always ensure you use the same size for both retrieval grid and gas VMRs!**

While the `gas_name` field accepts any string, we highly recommend to use ordinary molecule descriptors, such as `CO2`, `CH4`, and so on. Internally, there may be certain checks for these molecule strings. This is especially true for water vapor (`H2O` or `h2o`), which may be used to convert between specific humidity and a water vapor gas profile.

The `unit` field may be used to defile VMR profiles on *convenient* scales; e.g. methane profiles can be supplied in parts per billion (`u"ppb"`) while carbon dioxide profiles can be supplied in parts per million (`u"ppm"`). Generally, gas absorbers with different units can be used within the same atmosphere object, as the relevant calculations perform the unit conversions for each gas object separately.

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["gas_types.jl"]
Order = [:type]
```