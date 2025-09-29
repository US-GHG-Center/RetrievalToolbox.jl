# Miscellaneous functions

Functions here are usually general "helper" or "utility" type functions that a likely used throughout different parts of RetrievalToolbox.

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["misc.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```