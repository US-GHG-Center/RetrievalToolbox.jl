# Inversion functions

Here are listed all the generic inversion functions that act on any type of `AbstractSolver`. Ideally, users should use these functions to calculate or extract the various quantities of interest from the solver objects. In doing so, it is more likely that users can then more easily switch out solvers without having to adjust other parts of their algorithm.

## List of generic inversion functions

```@autodocs
Modules = [RetrievalToolbox]
Pages = ["inversion.jl"]
Order = [:function]
Filter = f -> !startswith(String(nameof(f)), "_") # Skip internal-use functions
```