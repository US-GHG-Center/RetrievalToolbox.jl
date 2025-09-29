# Atmosphere functions

Functions below are designed to operate on atmosphere object, and for the time being
mostly on `EarthAtmosphere` ones.

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["atmosphere.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```