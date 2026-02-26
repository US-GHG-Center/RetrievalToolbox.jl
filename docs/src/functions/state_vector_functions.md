# State Vector Functions

Below functions act on state vector elements (`AbstractStateVectorElement`), and (generally) nothing else. Their main purpose is mostly provide context-specific information about a given state vector element that may be needed to trigger other computations.

For example, the function `calculate_jacobian_before_isrf` is defined for all state vector elements (default implementation returns `true`) to indicate whether the Jacobians needs to be calculated before any instrument response functions are applied. The `ZeroLevelOffsetPolynomialSVE` type for example, represents an additive radiance contribution that can be added to the model radiance **after** the application instrument response function(s). While the result does not change, omitting the instrument function application for these state vector elements will reduce the overall computational effort.

Similarly, the function `is_aerosol_SVE` (default: `false`) tells the underlying radiative transfer code whether this particular state vector element is related to aerosols, and thus requires different handling in the allocation of linearized inputs to the RT solver.


## Aerosol Height
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_aerosol_height.jl"]
Order = [:function]
Filter = DocFilter
```

## Aerosol Optical Depth
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_aerosol_od.jl"]
Order = [:function]
Filter = DocFilter
```

## Aerosol Width
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_aerosol_width.jl"]
Order = [:function]
Filter = DocFilter
```

## BRDF Kernel Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_brdf_polynomial.jl"]
Order = [:function]
Filter = DocFilter
```

## Dispersion Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_dispersion.jl"]
Order = [:function]
Filter = DocFilter
```

## Gas VMR Profile
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_gas_profile.jl"]
Order = [:function]
Filter = DocFilter
```

## Gas Level Scaling Factor
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_gas_scale.jl"]
Order = [:function]
Filter = DocFilter
```

## ILS Stretch Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_ils_stretch.jl"]
Order = [:function]
Filter = DocFilter
```

## SIF Radiance
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_sif_radiance.jl"]
Order = [:function]
Filter = DocFilter
```

## Solar Scaler Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_solar_scale.jl"]
Order = [:function]
Filter = DocFilter
```

## Surface Albedo Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_surface_albedo.jl"]
Order = [:function]
Filter = DocFilter
```

## Surface Pressure
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_surface_pressure.jl"]
Order = [:function]
Filter = DocFilter
```

## Temperature Offset
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_temperature_offset.jl"]
Order = [:function]
Filter = DocFilter
```

## Zero Level Offset Polynomial
```@autodocs
Modules = [RetrievalToolbox]
Pages = ["state_vector_zero_level_offset.jl"]
Order = [:function]
Filter = DocFilter
```
